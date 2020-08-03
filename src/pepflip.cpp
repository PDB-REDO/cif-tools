/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 19 februari, 2018
*/

// test 3fvl

#include "pdb-redo.h"

#include <fcntl.h>

#include <iomanip>
#include <numeric>
#include <future>

#include <boost/program_options.hpp>


#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>

#include "cif++/Secondary.h"
#include "cif++/Statistics.h"
#include "cif++/CifUtils.hpp"

#include "minimizer.h"
#include "ramachandran.h"
#include "skiplist.h"

using namespace std;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;

using mmcif::Atom;
using mmcif::Point;
using mmcif::Structure;
using mmcif::Monomer;
using mmcif::BondMap;
using mmcif::Polymer;
using mmcif::StatsCollector;

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Coord_map;
using clipper::Coord_frac;

typedef mmcif::MapMaker<float> MapMaker;

// --------------------------------------------------------------------

size_t
	NTHREADS = boost::thread::hardware_concurrency();

float
	gMinAngle		= 90,
	gMinRSCC		= 0.6,
	gMaxRSCCDrop	= 0.95;

const double
	kMinRSCC = 0.6;

// --------------------------------------------------------------------

struct PFRes
{
	PFRes(int seqID, string compoundID, string entityID, string asymID, string authID, bool tryFlip)
		: mSeqID(seqID), mCompoundID(compoundID), mEntityID(entityID), mAsymID(asymID), mAuthID(authID), mTryFlip(tryFlip) {}
	
	int			mSeqID;
	string		mCompoundID, mEntityID, mAsymID, mAuthID;
	float		mTryAngle = 180;
	bool		mTrust = false;
	bool		mTryFlip;
	bool		mFlipped = false;
};

ostream& operator<<(ostream& os, const PFRes& r)
{
	os << r.mCompoundID << ' ' << r.mAsymID << r.mSeqID << " (" << r.mAuthID << ')';
	return os;
}

// --------------------------------------------------------------------

ostream& operator<<(ostream& os, const Atom& a);

ostream& operator<<(ostream& os, const RamachandranScore& s)
{
	switch (s)
	{
		case rsNotAllowed:		os << "N"; break;
		case rsAllowed:			os << "A"; break;
		case rsFavoured:		os << "F"; break;
	}

	return os;
}

// --------------------------------------------------------------------

class AtomLocationSaver
{
  public:

	AtomLocationSaver(mmcif::Residue& r)
		: mCommitted(false)
	{
		for (auto& a: r.atoms())
			mAtoms.emplace_back(make_pair(a, a.location()));
	}

	AtomLocationSaver(mmcif::Structure& s, const string& asymID, int seqFirst, int seqLast)
		: mCommitted(false)
	{
		for (auto& a: s.atoms())
		{
			if (a.labelAsymId() != asymID or a.labelSeqId() < seqFirst or a.labelSeqId() > seqLast)
				continue;
			
			mAtoms.emplace_back(make_pair(a, a.location()));
		}
	}
	
	~AtomLocationSaver()
	{
		if (not mCommitted)
		{
			for (auto& a: mAtoms)
				a.first.location(a.second);
		}
	}
	
	void commit() { mCommitted = true; }
	
	void storeAtomPositions(vector<pair<string,Point>>& atoms) const
	{
		atoms.reserve(mAtoms.size());
		for (auto& a: mAtoms)
			atoms.push_back(make_pair(a.first.id(), a.first.location()));
	}
	
  private:
	vector<pair<mmcif::Atom,mmcif::Point>> mAtoms;
	bool mCommitted;
};

// --------------------------------------------------------------------

void TrustBasedOnSS(Structure& structure, vector<PFRes>& residues)
{
	mmcif::DSSP dssp(structure);
	
	const set<mmcif::SecondaryStructureType> kTrustedDSSPTypes = {
			mmcif::ssAlphahelix, mmcif::ssHelix_3, mmcif::ssHelix_5, mmcif::ssStrand };
	const int kBufferSS = 1;
	
	enum { LOOP, SS_RUN } state = LOOP;
	mmcif::SecondaryStructureType cur = mmcif::ssLoop;
	
	size_t i = 0, b = 0;
	while (i < residues.size())
	{
		auto ss = dssp(residues[i].mAsymID, residues[i].mSeqID);
		++i;
		
		switch (state)
		{
			case LOOP:
				if (kTrustedDSSPTypes.count(ss))
				{
					state = SS_RUN;
					cur = ss;
					b = i - 1;
				}
				break;
			
			case SS_RUN:
				if (ss != cur)
				{
					// trust range, if long enough
					if (i - b > 2 * kBufferSS + 2)
					{
						for (auto j = b + kBufferSS; j < i - kBufferSS - 2; ++j)
						{
							if (cif::VERBOSE)
								cerr << residues[j] << "  - trusted based on secondary structure (" << char(cur) << ')' << endl;
							if (not dssp.isAlphaHelixEndBeforeStart(residues[j].mAsymID, residues[j].mSeqID))
								residues[j].mTrust = true;
						}

						if (cif::VERBOSE)
							cerr << residues[i - kBufferSS - 2] << "  - n-term of secondary structure (" << char(cur) << "). It will be considered for flips" << endl;
					}
					
					state = LOOP;
					--i;	// re-start
				}
				break;
		}
	}
}

class Rotation
{
  public:
	Rotation(Point ca1, Point ca2, float angle)
	{
		mT = ca1;
		
		angle = angle * mmcif::kPI / 180;
		
		auto s = sin(angle / 2);
		auto v = ca2 - ca1;
		v.normalize();
		
		mQ = mmcif::quaternion(cos(angle / 2), v.mX * s, v.mY * s, v.mZ * s);
	}

	Point operator()(const Point& p) const
	{
		Point result = p - mT;
		result.rotate(mQ);
		return result + mT;
	}

  private:
	Point mT;
	mmcif::quaternion mQ;
};

tuple<double,double,float,float> ScorePotentialFlips(MapMaker& mm, const Monomer& r1, const Monomer& r2, initializer_list<float> angles)
{
	double bestScoreB = 0, bestScoreD = 0;
	float bestAngleB = 0, bestAngleD = 0;
	
	auto& fb = static_cast<clipper::Xmap<float>&>(mm.fb());
	auto& fd = static_cast<clipper::Xmap<float>&>(mm.fd());
	auto& cell = mm.cell();
	
	for (auto angle: angles)
	{
		auto ca1 = r1.atomByID("CA");
		auto ca2 = r2.atomByID("CA");
		
		Rotation rt(ca1.location(), ca2.location(), angle);

		auto c = r1.atomByID("C");
		auto o = r1.atomByID("O");
		auto n = r2.atomByID("N");
//		auto h = r2.atomByID("H");

		double scoreB = 0, scoreD = 0;
		
		for (auto& a: { c, o, n })
		{
			auto l = rt(a.location());

			clipper::Coord_orth p = l;
			clipper::Coord_frac pf = p.coord_frac(cell);
			
			scoreB += static_cast<int>(a.type()) * a.occupancy() * fb.interp<clipper::Interp_cubic>(pf);
			scoreD += static_cast<int>(a.type()) * a.occupancy() * fd.interp<clipper::Interp_cubic>(pf);
		}

		if (bestScoreB < scoreB)
		{
			bestScoreB = scoreB;
			bestAngleB = angle;
		}

		if (bestScoreD < scoreD)
		{
			bestScoreD = scoreD;
			bestAngleD = angle;
		}
	}
	
	return make_tuple(bestScoreB, bestScoreD, bestAngleB, bestAngleD);
}

// --------------------------------------------------------------------

void CheckDensity(Structure& structure, vector<PFRes>& residues, MapMaker& mm)
{
	const double kMinRatioOrigFlip = 0.9;

	size_t potentialFlipCount = 0;
	
	auto& polymers = structure.polymers();
	
	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		auto& cur = residues[i];
		auto& next = residues[i + 1];

		if (next.mAsymID != cur.mAsymID)
		{
			cur.mTryFlip = false;
			continue;
		}

		string asymID = cur.mAsymID;
		auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
		assert(poly != polymers.end());
		
		auto& r1 = poly->getBySeqID(cur.mSeqID);
		auto& r2 = poly->getBySeqID(next.mSeqID);

		if (cur.mTrust or not cur.mTryFlip)
		{
			cout << cur << "  - residue excluded from peptide flips" << endl;
			cur.mTryFlip = false;
			continue;
		}

		cur.mTryFlip = false;
		
		if (not Monomer::areBonded(r1, r2))
		{
			cout << cur << "  - Distance between this CA and the next CA doesn't correspond with a peptide bond" << endl;
			continue;
		}
		
		try
		{
			double scoreBBefore, scoreBAfter, scoreDBefore, scoreDAfter;
			float angleBAfter, angleDAfter;

			tie(scoreBBefore, scoreDBefore, ignore, ignore) = ScorePotentialFlips(mm, r1, r2, { 330, 0, 30 });
			tie(scoreBAfter, scoreDAfter, angleBAfter, angleDAfter) = ScorePotentialFlips(mm, r1, r2, { 150, 180, 210 });
			
			if (scoreBAfter > kMinRatioOrigFlip * scoreBBefore)
			{
				cout << cur << "  - might need to be flipped based on the density fit, before: " << scoreBBefore << " after: " << scoreBAfter << endl;
				cur.mTryFlip = true;
				cur.mTryAngle = angleBAfter;
				++potentialFlipCount;
			}
			else if (scoreDAfter > scoreDBefore)
			{
				cout << cur << "  - might need to be flipped based on the difference density fit, before: " << scoreDBefore << " after: " << scoreDAfter << endl;
				cur.mTryFlip = true;
				cur.mTryAngle = angleDAfter;
				++potentialFlipCount;
			}
			else
				cout << cur << "  - no reason to flip this peptide based on the density" << endl;
		}
		catch (const exception& ex)
		{
			cerr << "Error processing " << cur << endl
				 << ex.what() << endl;	
		}
	}
}

// --------------------------------------------------------------------

//const regex kNotAHBondEnergyTypeRX("N(?:R(?:5|55|56|6|66))?");	// <-- N, NR5, NR55, NR56, NR6 or NR66
//
//int CalculateHBondsForAtoms(Structure& s, initializer_list<Atom> atoms)
//{
//	int result = 0;
//	
//	for (auto& a: atoms)
//	{
//		Polymer p(s, a.labelAsymId());
//		auto m = p.getBySeqID(a.labelSeqId());
//		
//		// all atoms near 3.5
//		for (auto& a2: dm.near(a, 3.5f))
//		{
//			// but only O or N
//			if (a2.type() != mmcif::O and a2.type() != mmcif::N)
//				continue;
//			
//			// And if N, not with energytype N, NR5, etc
//			if (a2.type() == mmcif::N)
//			{
//				auto& c = a2.comp();
//				const auto& ca = c.getAtomById(a2.labelAtomId());
//				if (regex_match(ca.typeEnergy, kNotAHBondEnergyTypeRX))
//					continue;
//			}
//			
//			// Not on same residue of course
//			if (a.labelSeqId() == a2.labelSeqId() and a.labelAsymId() == a2.labelAsymId())
//				continue;
//			
//			// if on same chain and both backbone, at least 3 residues apart
//			if (a.labelAsymId() == a2.labelAsymId() and a2.isBackBone())
//			{
//				auto m2 = p.getBySeqID(a2.labelSeqId());
//				if (p.Distance(m, m2) < 3)
//					continue;
//			}
//
//			if (cif::VERBOSE)
//				cerr << "HBond between " << a << " and " << a2 << endl;
//			
//			result += 1;
//		}
//	}
//	
//	return result;
//}

// --------------------------------------------------------------------

struct PepFlipScore
{
	size_t			fpix;
	string			id;
	float			angle;
	mmcif::Point	c;
	bool			refined;
	struct {
		double score;					// from minimizer
		double rsccsRes;
		float densityDifferenceO;
		RamachandranScore z1, z2;
//		int hbonds;
	} orig, flip;
	vector<pair<string,Point>> atoms;

	bool improved() const
	{
		return
				refined
			and
				abs(angle) > gMinAngle
			and
				flip.score < orig.score
			and
				flip.z1 >= rsAllowed and flip.z2 >= rsAllowed
			and
				(flip.z1 + flip.z2) >= (orig.z1 + orig.z2)
			and
				flip.rsccsRes > gMinRSCC
			and
				flip.rsccsRes >= gMaxRSCCDrop * orig.rsccsRes
			and	
				flip.densityDifferenceO > 0
			and
				flip.densityDifferenceO >= orig.densityDifferenceO
//			and
//				flip.hbonds >= orig.hbonds
			;
	}
	
	bool operator<(const PepFlipScore& rhs) const { return fpix < rhs.fpix; }
};

template<typename T>
inline auto best(const T& v, bool best, int prec)
{
	stringstream s;
	s << fixed << setprecision(prec) << setw(8) << v;
	string vs = s.str();
	return best ?
		coloured(vs.c_str(), cif::scCYAN, cif::scBLUE) :
		coloured(vs.c_str(), cif::scNONE, cif::scNONE, false);
}

template<>
inline auto best(const RamachandranScore& v, bool best, int prec)
{
	stringstream s;
	s << v;
	string vs = s.str();
	return best ?
		coloured(vs.c_str(), cif::scCYAN, cif::scBLUE) :
		coloured(vs.c_str(), cif::scNONE, cif::scNONE, false);
}

template<>
inline auto best(const int& v, bool best, int prec)
{
	string s = to_string(v);
	return best ?
		coloured(s.c_str(), cif::scCYAN, cif::scBLUE) :
		coloured(s.c_str(), cif::scNONE, cif::scNONE, false);
}

struct headers1 {};
struct headers2 {};

ostream& operator<<(ostream& os, const headers1&)
{
	os << "Angle    "
	   << "Score             "
	   << "RSCCS residues    "
	   << "Diff.Density O    "
	   << "RamaZ  "
	   ;
	return os;
}

ostream& operator<<(ostream& os, const headers2)
{
	os << "         "
	   << "Original Flipped  "
	   << "Original Flipped  "
	   << "Original Flipped  "
	   << "O F O F "
	   ;
	return os;
}

ostream& operator<<(ostream& os, const PepFlipScore& rhs)
{
	os << best(rhs.angle, abs(rhs.angle) > 90.0f, 1) << ' '

	   << best(rhs.orig.score, false, 1) << ' '
	   << best(rhs.flip.score, rhs.flip.score < rhs.orig.score, 1) << ' '

	   << best(rhs.orig.rsccsRes, false, 3) << ' '
	   << best(rhs.flip.rsccsRes, rhs.flip.rsccsRes > rhs.orig.rsccsRes, 3) << ' '

	   << best(rhs.orig.densityDifferenceO, false, 3) << ' '
	   << best(rhs.flip.densityDifferenceO, rhs.flip.densityDifferenceO > 0 and rhs.flip.densityDifferenceO > rhs.orig.densityDifferenceO, 3) << ' '

	   << best(rhs.orig.z1, false, 0) << ' '
	   << best(rhs.flip.z1, rhs.flip.z1 > rhs.orig.z1, 0) << ' '

	   << best(rhs.orig.z2, false, 0) << ' '
	   << best(rhs.flip.z2, rhs.flip.z2 > rhs.orig.z2, 0) << ' '

//	   << best(rhs.orig.hbonds, false, 0) << ' '
//	   << best(rhs.flip.hbonds, rhs.flip.hbonds > rhs.orig.hbonds, 1)
	   
	   ;

	return os;
}

PepFlipScore FlipPeptide(Structure& structure, const Polymer& poly,
	const vector<PFRes>& residues, size_t i, StatsCollector& sc, MapMaker& mm,
	mmcif::BondMap& bm, const string& algorithm, float mapWeight, float plane5AtomsESD,
	bool testMode)
{
	auto& r1 = poly.getBySeqID(residues[i].mSeqID);
	auto& r2 = poly.getBySeqID(residues[i + 1].mSeqID);
	
	const clipper::Xmap<float>& fd = mm.fd();
	
	auto getDf = [&fd](const Point& p)
	{
		Coord_orth po = p;
		Coord_frac pf = po.coord_frac(fd.cell());
		return fd.interp<clipper::Interp_cubic>(pf);
	};
	
	AtomLocationSaver s(structure, r1.asymID(), r1.seqID() - 1, r2.seqID() + 1);

	PepFlipScore result = { i, boost::lexical_cast<string>(residues[i]) };

	auto prePro = [&](size_t i)
	{
		return (i + 1 < residues.size() and residues[i + 1].mCompoundID == "PRO" and residues[i].mCompoundID != "PRO");
	};

	auto ca1 = r1.atomByID("CA");
	auto ca2 = r2.atomByID("CA");
	
	auto c = r1.atomByID("C");
	auto o = r1.atomByID("O");
	auto n = r2.atomByID("N");
//	auto h = r2.atomByID("H");

	result.c = c.location();

	unique_ptr<Minimizer> origMinimizer(Minimizer::create(algorithm, poly, r1.seqID() - 1, r2.seqID() + 1, bm, mm.fb(), mapWeight, plane5AtomsESD));
	result.orig.score = origMinimizer->refine(false);

	result.orig.rsccsRes = sc.collect({ &r1, &r2 }).RSCCS;
	result.orig.densityDifferenceO = getDf(o.location());
	result.orig.z1 = calculateRamachandranScore(r1.compoundID(), prePro(i), r1.phi(), r1.psi());
	result.orig.z2 = calculateRamachandranScore(r2.compoundID(), prePro(i + 1), r2.phi(), r2.psi());
//	result.orig.hbonds = CalculateHBondsForAtoms(structure, { o, n }, dm);

//	auto cl = c.location();
	auto ol = o.location();
//	auto nl = n.location();

	auto angle = residues[i].mTryAngle;
	Rotation rt(ca1.location(), ca2.location(), angle);

	c.location(rt(c.location()));
	o.location(rt(o.location()));
	n.location(rt(n.location()));

	if (testMode)
		structure.getFile().save("unrefined-" + r1.labelID() + ".pdb");

	unique_ptr<Minimizer> flipMinimizer(Minimizer::create(algorithm, poly, r1.seqID() - 1, r2.seqID() + 1, bm, mm.fb(), mapWeight, plane5AtomsESD));

	double unrefinedScore = flipMinimizer->score();
	result.flip.score = flipMinimizer->refine(true);
	result.refined = result.flip.score < unrefinedScore;

	if (testMode)
		structure.getFile().save("refined-" + r1.labelID() + ".pdb");

//	auto rcl = c.location();
	auto rol = o.location();
//	auto rnl = n.location();

	result.angle = DihedralAngle(ol, ca1.location(), ca2.location(), rol);

	result.flip.rsccsRes = sc.collect({ &r1, &r2 }).RSCCS;
	result.flip.densityDifferenceO = getDf(o.location());
	result.flip.z1 = calculateRamachandranScore(r1.compoundID(), prePro(i), r1.phi(), r1.psi());
	result.flip.z2 = calculateRamachandranScore(r2.compoundID(), prePro(i + 1), r2.phi(), r2.psi());
//	result.flip.hbonds = CalculateHBondsForAtoms(structure, { o, n });

//cout << residues[i] << endl;
	if (result.improved())
		result.atoms = flipMinimizer->getAtoms();

	return result;
}

vector<size_t> CombineTrustBasedOnNCS(const list<Polymer>& polymers, vector<PFRes>& residues)
{
	typedef vector<PFRes>::iterator PFResIter;
	
	set<string> entities;
	for (auto& p: polymers)
		entities.insert(p.entityID());
	
	vector<size_t> result;
	
	for (const string& entityID: entities)
	{
		vector<pair<PFResIter,PFResIter>> pri;
		
		for (auto ri = residues.begin(); ri != residues.end(); ++ri)
		{
			if (ri->mEntityID == entityID)
			{
				auto bri = ri, eri = ri + 1;
				
				while (eri != residues.end() and eri->mEntityID == entityID and eri->mAsymID == bri->mAsymID)
					++eri;
				
				pri.emplace_back(bri, eri);
				ri = eri - 1;
			}
		}

		int seqID = 0;
		vector<PFResIter> tri;

		while (pri.size() > 1)
		{
			tri.clear();
			for (auto& ri: pri)
			{
				if (ri.first != ri.second and ri.first->mSeqID == seqID)
				{
					tri.push_back(ri.first);
					++ri.first;
				}
			}
			
			if (tri.size() < 2)
			{
				++seqID;
				continue;
			}
			
			size_t flipCount = accumulate(tri.begin(), tri.end(), 0, [](int sum, auto ri) { return sum + (ri->mFlipped ? 1 : 0); });
			if (flipCount < tri.size() and flipCount > 0)
			{
				for (auto ri: tri)
				{
					if (ri->mTryFlip)
						continue;
					
					result.push_back(ri - residues.begin());
				}
			}
			
			pri.erase(remove_if(pri.begin(), pri.end(), [](auto i) { return i.first == i.second; }), pri.end());
		}
	}
	
	return result;
}

void JoinFlips(Structure& structure, const Polymer& poly, MapMaker& mm, mmcif::BondMap& bm,
	const string& algorithm, float mapWeight, float plane5AtomsESD,
	const vector<PFRes>& residues, vector<PepFlipScore>& flipped, size_t f1, size_t f2)
{
	size_t ri1 = flipped[f1].fpix;
	size_t ri2 = flipped[f2].fpix;
	
	int s1 = residues[ri1].mSeqID - 1;
	int s2 = residues[ri2].mSeqID + 1;
	
	cout << "Joining flips from " << residues[ri1] << " to " << residues[ri2] << endl;
	
	// Eerste poging, gewoon alles flippen en dan verfijnen...
	
	for (size_t fi = f1; fi <= f2; ++fi)
	{
		size_t ri = flipped[fi].fpix;
		
		auto& r1 = poly.getBySeqID(residues[ri].mSeqID);
		auto& r2 = poly.getBySeqID(residues[ri + 1].mSeqID);

		auto ca1 = r1.atomByID("CA");
		auto ca2 = r2.atomByID("CA");
		
		auto c = r1.atomByID("C");
		auto o = r1.atomByID("O");
		auto n = r2.atomByID("N");
		
		auto angle = residues[ri].mTryAngle;
		Rotation rt(ca1.location(), ca2.location(), angle);
	
		c.location(rt(c.location()));
		o.location(rt(o.location()));
		n.location(rt(n.location()));
		
		flipped[fi].atoms.clear();
	}
	
	unique_ptr<Minimizer> flipMinimizer(Minimizer::create(algorithm, poly,
		s1, s2, bm, mm.fb(), mapWeight, plane5AtomsESD));

	double unrefinedScore = flipMinimizer->score();
	double refinedScore = flipMinimizer->refine(true);

	if (refinedScore > unrefinedScore)
	{
		cerr << "Oeps, deze verfijnde niet" << endl;
		
		for (size_t fi = f1; fi <= f2; ++fi)
			flipped[fi].refined = false;
	}
}

void FlipPeptides(Structure& structure, const string& asymID,
	int resFirst, int resLast, const SkipList& skip,
	MapMaker& mm, bool trustDSSP,
	const string& algorithm, float mapWeight, float plane5AtomsESD,
	ofstream& cootScript, bool testMode)
{
	StatsCollector sc(mm, structure, false /*electronScattering*/);

	auto& polymers = structure.polymers();
	
	vector<PFRes> residues;
	for (auto& p: polymers)
	{
		for (auto& m: p)
		{
			bool tryFlip = asymID.empty() or (m.asymID() == asymID and m.seqID() >= resFirst and m.seqID() <= resLast);
			
			if (tryFlip)
				tryFlip = find_if(skip.begin(), skip.end(), [&](auto s) { return s.asymID == m.asymID() and s.seqID == m.seqID(); }) == skip.end();
			
			residues.emplace_back(m.seqID(), m.compoundID(), p.entityID(), m.asymID(), m.authID(), tryFlip);
		}
	}
	
	if (trustDSSP)
	{
		if (cif::VERBOSE)
			cerr << endl
				 << "Calculating DSSP" << endl;
	
		TrustBasedOnSS(structure, residues);

		if (cif::VERBOSE)
			cerr << "DSSP done" << endl;
	}

	if (cif::VERBOSE)
		cerr << endl
			 << "Initial check" << endl;

	auto df = async(launch::async, bind(CheckDensity, ref(structure), ref(residues), ref(mm)));

	mmcif::BondMap bm(structure);
	
	df.wait();

	cout << endl
	     << "Potential number of flips: " << accumulate(residues.begin(), residues.end(), 0, [](int sum, auto& r) { return sum + (r.mTryFlip ? 1 : 0); }) << endl
	     << endl;

	vector<PepFlipScore> flipped;
	flipped.reserve(residues.size());
	
	unique_ptr<cif::Progress> p;
	if (not cif::VERBOSE)
		p.reset(new cif::Progress(residues.size(), "Flipping, 1st round"));
	
	if (NTHREADS == 1)
	{
		for (size_t i = 0; i + 1 < residues.size(); ++i)
		{
			if (p)
				p->progress(i);
	
			if (residues[i].mTrust or not residues[i].mTryFlip)
				continue;
	
			string asymID = residues[i].mAsymID;
			auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
			assert(poly != polymers.end());
	
			try
			{
				flipped.emplace_back(
					FlipPeptide(structure, *poly, residues, i, sc, mm, bm, algorithm, mapWeight, plane5AtomsESD, testMode));
				residues[i].mFlipped = flipped.back().improved();
			}
			catch (const exception& ex)
			{
				cerr << "Error processing " << residues[i] << ": " << ex.what() << endl;
				continue;
			}
		}
	}
	else
	{
		boost::thread_group t;
		atomic<size_t> next(0);
		mutex m;
		
		for (size_t ti = 0; ti < NTHREADS; ++ti)
			t.create_thread([&]()
			{
				Structure s(structure);
				StatsCollector sc(mm, s, false /*electronScattering*/);
				auto& polymers = s.polymers();

				for (;;)
				{
					size_t i = next++;
					
					if (i + 1 >= residues.size())
						break;
					
					if (p)
						p->progress(i);
					
					if (residues[i].mTrust or not residues[i].mTryFlip)
						continue;
			
					string asymID = residues[i].mAsymID;
					auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
					assert(poly != polymers.end());
			
					try
					{
						auto flipResult = FlipPeptide(s, *poly, residues, i, sc, mm, bm, algorithm, mapWeight, plane5AtomsESD, testMode);
						residues[i].mFlipped = flipResult.improved();
						
						unique_lock<mutex> lock(m);
						flipped.emplace_back(move(flipResult));
					}
					catch (const exception& ex)
					{
						cerr << "Error processing " << residues[i] << ": " << ex.what() << endl;
						continue;
					}
				}
			});
	
		t.join_all();
		
		sort(flipped.begin(), flipped.end());
		
	}

	if (p)
		p->progress(residues.size());

	cout << endl;

	// second round, flip any residue whose NCS copy was flipped
	if (asymID.empty())
	{
		auto residuesRound2 = CombineTrustBasedOnNCS(polymers, residues);
		if (not residuesRound2.empty())
		{
			cout << endl
				 << "Potential number of flips in second round (NCS copies): " << residuesRound2.size() << endl
				 << endl;
			
			unique_ptr<cif::Progress> p;
			if (not cif::VERBOSE)
				p.reset(new cif::Progress(residuesRound2.size(), "Flipping, 2nd round"));
			
			for (size_t i: residuesRound2)
			{
				string asymID = residues[i].mAsymID;
				auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
				assert(poly != polymers.end());
		
				try
				{
					flipped.emplace_back(
						FlipPeptide(structure, *poly, residues, i, sc, mm, bm, algorithm, mapWeight, plane5AtomsESD, testMode));
					residues[i].mFlipped = flipped.back().improved();
				}
				catch (const exception& ex)
				{
					cerr << "Error processing " << residues[i] << ": " << ex.what() << endl;
					continue;
				}
	
				if (p)
					p->consumed(1);
			}
		}
	}
	
	// third step, some of the flips might overlap, process them separately
	for (size_t i = 0; i + 1 < flipped.size(); ++i)
	{
		if (not flipped[i].improved())
			continue;
		
		string asymID = residues[flipped[i].fpix].mAsymID;
		int firstSeq = residues[flipped[i].fpix].mSeqID - 1;
		int lastSeq = firstSeq + 3;
		
		size_t j = i;
		while (j + 1 < flipped.size() and flipped[j + 1].improved() and lastSeq + 1 >= residues[flipped[j + 1].fpix].mSeqID)
		{
			lastSeq = residues[flipped[j].fpix].mSeqID + 2;
			++j;
		}
		
		if (j > i)
		{
			auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
			JoinFlips(structure, *poly, mm, bm, algorithm, mapWeight, plane5AtomsESD, residues, flipped, i, j);
		}
	}

	size_t idLen = 12;
	for (auto& s: flipped)
		if (idLen < s.id.length())
			idLen = s.id.length();
	
	cout << endl
		 << string(idLen + 8, ' ') << headers1() << endl
		 << string(idLen + 8, ' ') << headers2() << endl
		 << string(cif::get_terminal_width(), '-') << endl;
	
	vector<string> flippedIDs;
	for (auto& score: flipped)
	{
		string id = score.id + string(idLen - score.id.length(), ' ');

		auto red = [](const char* s){ return cif::coloured(s, cif::scWHITE, cif::scRED); };
		
		if (not score.improved())
		{
			cout << id << " " << "noflip" << " " << score << endl;
			continue;
		}

		cout << id << " " << red("flip") << "   " << score << endl;

		if (cootScript.is_open())
			cootScript << "(list \"Flipped " << score.id << "\" "
					   << score.c.mX << ' ' << score.c.mY << ' ' << score.c.mZ << " )" << endl;
		
		for (auto a: score.atoms)
			structure.getAtomById(a.first).location(a.second);
		
		flippedIDs.push_back(score.id);
	}

	cout << endl
		 << string(cif::get_terminal_width(), '-') << endl
		 << "Summary: Flipped " << flippedIDs.size() << " peptides" << endl
		 << ba::join(flippedIDs, ", ") << endl;
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("hklin",				po::value<string>(),	"reflections file")
		("xyzin",				po::value<string>(),	"coordinates file")
		("xyzout,o",			po::value<string>(),	"output coordinates to this file")
		
		("log",					po::value<string>(),	"Write verbose output to log file")
		
		("asym-id",				po::value<string>(),	"Asymetric unit to pepflip")
		("algorithm",			po::value<string>(),	"Refinement algorithm to use, can be sa or gsl")
		("res-first",			po::value<int>(),		"Sequence number for first residue to pepflip, default = 1")
		("res-last",			po::value<int>(),		"Sequence number for last residue to pepflip, default is last in sequence")
		("dict",				po::value<vector<string>>(),
														"Dictionary file containing restraints for residues in this specific target, can be specified multiple times.")

		("no-dssp",										"By default, residues in helix and strand are trusted, use this switch to turn this off")

		("skip",				po::value<string>(),	"File containing the skip list: the residues that should not be flipped")

		("map-weight",			po::value<float>(),		"Map weight in minimisation, default is 60")

		("minimal-angle",		po::value<float>(),		"Minimal angle for flip in degrees")
		("minimal-rscc",		po::value<float>(),		"Minimal RSCC score for the two residues")
		("max-rscc-drop",		po::value<float>(),		"The flipped rscc should be at least this factor times the orignal rscc")

		("plane-5-atoms-esd",	po::value<float>(),		"ESD for the atoms in the 5 atom peptide bond plane, default is 0.11")
		("coot-script",									"Write a Coot script for the changed peptides")

		("help,h",										"Display help message")
		("version",										"Print version")
		
		("sampling-rate",		po::value<float>(),		"Sampling rate")
		("max-threads",			po::value<uint16_t>(),	"Max number of threads")
		
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)")
		("test",										"Do a test run, flipping all potential residues and writing intermediate files for further debugging")
		;

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("hklin", 1);
	p.add("xyzin", 1);
	p.add("xyzout", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "pepflip.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "pepflip.conf";
	
	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}
	
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		cerr << visible_options << endl;
		exit(0);
	}
	
	if (vm.count("xyzin") == 0 or vm.count("hklin") == 0)
	{
		cerr << "Input files not specified" << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	if (vm.count("log"))
	{
		string logFile = vm["log"].as<string>();
		
	    // open the log file
	    int fd = open(logFile.c_str(), O_CREAT|O_APPEND|O_RDWR, 0644);
	    if (fd < 0)
	        throw runtime_error("Opening log file " + logFile + " failed: " + strerror(errno));
	
	    // redirect stdout and stderr to the log file
	    dup2(fd, STDOUT_FILENO);
	    dup2(fd, STDERR_FILENO);
	    close(fd);
	}

	if (vm.count("max-threads"))
		NTHREADS = vm["max-threads"].as<uint16_t>();
	if (NTHREADS > boost::thread::hardware_concurrency())
		NTHREADS = boost::thread::hardware_concurrency();

	if (cif::VERBOSE or NTHREADS < 1)
		NTHREADS = 1;

	if (vm.count("dict"))
	{
		for (auto dict: vm["dict"].as<vector<string>>())
			mmcif::CompoundFactory::instance().pushDictionary(dict);
	}

	mmcif::File f(vm["xyzin"].as<string>());
	Structure structure(f);
	
	fs::path mtzFile = vm["hklin"].as<string>();

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

	MapMaker mm;
	mm.loadMTZ(mtzFile, samplingRate);

	string asymID;
	if (vm.count("asym-id"))
		asymID = vm["asym-id"].as<string>();
		
	int resFirst = 1;
	if (vm.count("res-first"))
		resFirst = vm["res-first"].as<int>();
	
	int resLast = numeric_limits<int>::max();
	if (vm.count("res-last"))
		resLast = vm["res-last"].as<int>();

	string algorithm = "gsl";
	if (vm.count("algorithm"))
		algorithm = vm["algorithm"].as<string>();

	if (vm.count("minimal-angle"))	gMinAngle = vm["minimal-angle"].as<float>();
	if (vm.count("minimal-rscc"))	gMinRSCC = vm["minimal-rscc"].as<float>();
	if (vm.count("max-rscc-drop"))	gMaxRSCCDrop = vm["max-rscc-drop"].as<float>();


	bool trustDSSP = vm.count("no-dssp") == 0;
	
	float mapWeight = 60;
	if (vm.count("map-weight"))
		mapWeight = vm["map-weight"].as<float>();
	
	float plane5AtomsESD = 0.11;
	if (vm.count("plane-5-atoms-esd"))
		plane5AtomsESD = vm["plane-5-atoms-esd"].as<float>();
	
	fs::path xyzout;
	if (vm.count("xyzout"))
		xyzout = vm["xyzout"].as<string>();
	else
	{
		xyzout = vm["xyzin"].as<string>();
		if (xyzout.extension() == ".gz")
			xyzout = xyzout.stem();
		xyzout = xyzout.parent_path() / (xyzout.filename().stem().string() + "-flipped.cif");
	}
	
	xyzout = fs::system_complete(xyzout);
	
	fs::ofstream cootScript;
	if (vm.count("coot-script"))
	{
		fs::path cs = xyzout.extension() == ".gz" ?
			xyzout.parent_path() / (xyzout.filename().stem().stem().string() + ".scm") :
			xyzout.parent_path() / (xyzout.filename().stem().string() + ".scm");
		cootScript.open(cs);
		
		if (not cootScript.is_open())	
			throw runtime_error("Failed to open coot script " + cs.string());
		
		cootScript << "(interesting-things-gui \"Pepflips\"" << endl
				   << "(list" << endl;
	}
	
	bool testMode = false;
	if (vm.count("test"))
	{
		testMode = true;
		
		fs::path testOutputDir = fs::current_path() / (xyzout.filename().stem().string() + "-test");
		
		if (not fs::exists(testOutputDir))
			fs::create_directory(testOutputDir);
		
		fs::current_path(testOutputDir);
		
		NTHREADS = 1;
	}
	
	SkipList skip;
	if (vm.count("skip"))
		skip = readSkipList(vm["skip"].as<string>(), structure);
	
	FlipPeptides(structure, asymID,
		resFirst, resLast, skip, mm, trustDSSP, algorithm, mapWeight, plane5AtomsESD, cootScript, testMode);
	
	if (cootScript.is_open())
	{
		cootScript << "))" << endl;
		cootScript.close();
	}

	f.save(xyzout);
	
	return 0;
}
