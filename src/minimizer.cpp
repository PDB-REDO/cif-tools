/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 22 mei, 2018
*/

#include "pdb-redo.h"

#include <future>
#include <regex>
#include <iomanip>
#include <random>

#include <boost/format.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/thread.hpp>

#include "minimizer.h"

using namespace std;
using namespace mmcif;
namespace fs = boost::filesystem;

// --------------------------------------------------------------------

const uint32_t kRefSentinel = numeric_limits<uint32_t>::max();

const double
	kNonBondedContactDistanceSq = 11.0 * 11.0,
	kMaxPeptideBondLengthSq = 3.5 * 3.5;

// --------------------------------------------------------------------

ostream& operator<<(ostream& os, const Atom& a)
{
	os << a.labelAsymId() << ':' << a.labelSeqId() << '/' << a.labelAtomId();
	
	return os;
}

// --------------------------------------------------------------------

struct lessAtom
{
	bool operator()(const Atom& a, const Atom& b) const { return a.id().compare(b.id()) < 0; }
};

typedef set<Atom, lessAtom> AtomSet;

// --------------------------------------------------------------------

DPoint AtomLocationProvider::operator[](AtomRef atomID) const
{
	if (atomID >= mAtoms.size())
		throw range_error("Unknown atom " + to_string(atomID));
	return mAtoms[atomID].location();
}

string AtomLocationProvider::atom(AtomRef atomID) const
{
	if (atomID >= mAtoms.size())
		throw range_error("Unknown atom " + to_string(atomID));
	auto& a = mAtoms[atomID];
	return to_string(a.labelSeqId()) + ' ' + a.labelAtomId();
}

// --------------------------------------------------------------------

Minimizer::Minimizer(const mmcif::Polymer& poly, int first, int last,
	const mmcif::BondMap& bonds, const XMap& xMap, float mapWeight, float plane5AtomsESD)
{
	const Monomer* prev = nullptr;				// used to link residues
	
	for (auto& r: poly)
	{
		if (r.seqID() < first)
		{
			prev = &r;
			continue;
		}

		if (r.seqID() <= last)
		{
			for (auto a: r.atoms())
			{
				(void)ref(a);
				mAtoms.push_back(a);
			}
		}
	
		if (prev != nullptr)
			try
			{
				Atom c = prev->atomByID("C"), n = r.atomByID("N");
	
				if (DistanceSquared(c, n) < kMaxPeptideBondLengthSq)
				{
					bool trans = not Monomer::isCis(*prev, r);
					
					if (trans)
						addLinkRestraints(*prev, r, r.compoundID() == "PRO" ? "PTRANS" : "TRANS");
					else
						addLinkRestraints(*prev, r, r.compoundID() == "PRO" ? "PCIS" : "CIS");
					
					if (trans)
					{
						Atom ca1 = prev->atomByID("CA");
						Atom ca2 = r.atomByID("CA");
						
						mTransPeptideRestraints.emplace_back(TransPeptideRestraint{
							ref(ca1),
							ref(c),
							ref(n),
							ref(ca2)
						});
					}
		
					// add planar restraints
					vector<AtomRef> atoms = {
						ref(prev->atomByID("CA")),
						ref(c),
						ref(prev->atomByID("O")),
						ref(n),
						ref(r.atomByID("CA"))
					};
	
					mPlanarityRestraints.emplace_back(PlanarityRestraint{move(atoms), plane5AtomsESD});
				}
			}
			catch (const exception& ex)
			{
				if (cif::VERBOSE)
					cerr << "While processing plane-5-atoms restraints: " << ex.what() << endl;
	//			continue;
			}
		
		if (r.seqID() > last)
			break;
	
		prev = &r;
				
		auto& compound = r.compound();

		for (auto& b: compound.bonds())
		{
			try
			{
				if (compound.getAtomById(b.atomID[0]).typeSymbol == H or
					compound.getAtomById(b.atomID[1]).typeSymbol == H)
				{
					continue;
				}
				
				Atom a1 = r.atomByID(b.atomID[0]);
				Atom a2 = r.atomByID(b.atomID[1]);

				mBondRestraints.emplace_back(ref(a1), ref(a2), b.distance, b.esd);
			}
			catch (const exception& ex)
			{
				if (cif::VERBOSE > 1)
					cerr << "While processing bond restraints: " << ex.what() << endl;
				continue;
			}
		}

		for (auto& a: compound.angles())
		{
			try
			{
				if (compound.getAtomById(a.atomID[0]).typeSymbol == H or
					compound.getAtomById(a.atomID[1]).typeSymbol == H or
					compound.getAtomById(a.atomID[2]).typeSymbol == H)
				{
					continue;
				}
				
				Atom a1 = r.atomByID(a.atomID[0]);
				Atom a2 = r.atomByID(a.atomID[1]);
				Atom a3 = r.atomByID(a.atomID[2]);

				mAngleRestraints.emplace_back(
					ref(a1), ref(a2), ref(a3), a.angle, a.esd);
			}
			catch (const exception& ex)
			{
				if (cif::VERBOSE > 1)
					cerr << "While processing angle restraints: " << ex.what() << endl;
				continue;
			}
		}
		
//		for (auto& a: compound.torsions())
//		{
//			if (a.esd == 0)
//				continue;
//		
//			try
//			{
//				if (compound.getAtomById(a.atomID[0]).typeSymbol == H or
//					compound.getAtomById(a.atomID[1]).typeSymbol == H or
//					compound.getAtomById(a.atomID[2]).typeSymbol == H or
//					compound.getAtomById(a.atomID[3]).typeSymbol == H)
//				{
//					continue;
//				}
//				
//				Atom a1 = r.atomByID(a.atomID[0]);
//				Atom a2 = r.atomByID(a.atomID[1]);
//				Atom a3 = r.atomByID(a.atomID[2]);
//				Atom a4 = r.atomByID(a.atomID[3]);
//
//				mTorsionRestraints.emplace_back(a1, a2, a3, a4, a.angle, a.esd, a.period);
//			}
//			catch (const exception& ex)
//			{
//				if (cif::VERBOSE > 1)
//					cerr << "While processing torsion restraints: " << ex.what() << endl;
//				continue;
//			}
//		}
		
		for (auto& cv: compound.chiralCentres())
		{
			try
			{
				if (compound.getAtomById(cv.atomID[0]).typeSymbol == H or
					compound.getAtomById(cv.atomID[1]).typeSymbol == H or
					compound.getAtomById(cv.atomID[2]).typeSymbol == H)
				{
					continue;
				}
				
				Atom cc = r.atomByID(cv.atomIDCentre);
				Atom a1 = r.atomByID(cv.atomID[0]);
				Atom a2 = r.atomByID(cv.atomID[1]);
				Atom a3 = r.atomByID(cv.atomID[2]);

				auto volume = compound.chiralVolume(cv.id);

				mChiralVolumeRestraints.emplace_back(ref(cc), ref(a1),
					ref(a2), ref(a3), volume * 6);
			}
			catch (const exception& ex)
			{
				if (cif::VERBOSE > 1)
					cerr << "While processing chiral volume restraints: " << ex.what() << endl;
				continue;
			}
		}
		
		for (auto& p: compound.planes())
		{
			try
			{
				vector<AtomRef> atoms;
				
				for (auto a: p.atomID)
				{
					if (compound.getAtomById(a).typeSymbol == H)
						continue;
					
					atoms.push_back(ref(r.atomByID(a)));
				}

				if (atoms.size() > 3)
					mPlanarityRestraints.emplace_back(PlanarityRestraint{move(atoms), p.esd});
			}
			catch (const exception& ex)
			{
				if (cif::VERBOSE > 1)
					cerr << "While processing planarity restraints: " << ex.what() << endl;
				continue;
			}
		}
	}

	// Add link BondRestraints
//	for (auto& a1: mAtoms)
//	{
//		// TODO: implement based on struct_conn information and radii from ener_lib?
//	}
	
	if (mAtoms.empty())
		throw runtime_error("No atoms to refine");

	fs::path enerLibFilePath(getenv("CLIBD_MON"));
	enerLibFilePath /= "ener_lib.cif";

	cif::File enerLibFile(enerLibFilePath);
	auto& db = enerLibFile.firstDatablock();
	auto& libAtom = db["lib_atom"];

	const regex donorRx("B|D|H"), acceptorRx("B|A|H");

	set<tuple<AtomRef,AtomRef>> nbc;

	// now add the non-bonded restraints
	for (auto& a1: mAtoms)
	{
		AtomRef ra1 = ref(a1);
		
		for (auto& a2: poly.structure()->atoms())
		{
			if (a1 == a2)
				continue;

			if (bonds(a1, a2))
				continue;

			if (DistanceSquared(a1, a2) > kNonBondedContactDistanceSq)
				continue;

			AtomRef ra2 = ref(a2);

			if (nbc.count(make_tuple(ra1, ra2)))
				continue;
			
			if (find_if(mAngleRestraints.begin(), mAngleRestraints.end(),
				[&](auto& ar) { return (ar.mA == ra1 and ar.mC == ra2) or (ar.mA == ra2 and ar.mC == ra1); }) != mAngleRestraints.end())
				continue;

			if ((a1.labelCompId() == "PRO" or a1.labelCompId() == "HYP") and
				a1.labelSeqId() == a2.labelSeqId() + 1 and
				a1.labelAtomId() == "CD")
			{
				continue;
			}

			if ((a2.labelCompId() == "PRO" or a2.labelCompId() == "HYP") and
				a2.labelSeqId() == a1.labelSeqId() + 1 and
				a2.labelAtomId() == "CD")
			{
				continue;
			}
			
			if ((a1.labelCompId() == "ASN" or a2.labelCompId() == "NAG") and
				a1.labelAtomId() == "OD1" and a2.labelAtomId() == "C1")
			{
				continue;
			}

			if ((a1.labelCompId() == "NAG" or a2.labelCompId() == "ASN") and
				a1.labelAtomId() == "C1" and a2.labelAtomId() == "OD1")
			{
				continue;
			}

			double minDist = 2.8;
			
			if (bonds.is1_4(a1, a2))
			{
				if (cif::VERBOSE > 1)
					cerr << "1_4 for " << a1 << " and " << a2 << endl;
				minDist = 2.64;
			}
			else if ((a1.labelSeqId() + 1 == a2.labelSeqId() and a1.labelAtomId() == "O" and a2.labelAtomId() == "C") or
					 (a2.labelSeqId() + 1 == a1.labelSeqId() and a2.labelAtomId() == "O" and a1.labelAtomId() == "C"))
			{
				minDist = 2.84;
			}
			else
			{
				try
				{
					string et1 = a1.energyType();
					string et2 = a2.energyType();
					
					if (not (et1.empty() or et2.empty()))
					{
						auto r1 = libAtom.find(cif::Key("type") == et1);
						auto r2 = libAtom.find(cif::Key("type") == et2);
						
						if (not (r1.empty() or r2.empty()))
						{
							if (AtomTypeTraits(a1.type()).isMetal())
								minDist = r1.front()["ion_radius"].as<float>();
							else
								minDist = r1.front()["vdw_radius"].as<float>();
							
							if (AtomTypeTraits(a2.type()).isMetal())
								minDist += r2.front()["ion_radius"].as<float>();
							else
								minDist += r2.front()["vdw_radius"].as<float>();

							// OK, now that we're here, see if the atoms are in the same residue...
							
							if (a1.labelAsymId() == a2.labelAsymId() and a1.labelSeqId() == a2.labelSeqId())
								minDist *= 0.84;
							
							string hbType1 = r1.front()["hb_type"].as<string>(),
								   hbType2 = r2.front()["hb_type"].as<string>();
							
							if (regex_match(hbType1, donorRx) and regex_match(hbType2, acceptorRx))
							{
								minDist -= 0.5;
								if (hbType1 == "H")
									minDist -= 0.3;
							}

							if (regex_match(hbType2, donorRx) and regex_match(hbType1, acceptorRx))
							{
								minDist -= 0.5;
								if (hbType2 == "H")
									minDist -= 0.3;
							}
						}
					}
					
					// so-called strange exceptions in coot code
					
					if (find(mAtoms.begin(), mAtoms.end(), a2) == mAtoms.end())
					{
						switch (abs(a1.labelSeqId() - a2.labelSeqId()))
						{
							case 1:
								if ((a1.labelAtomId() == "O" and a2.labelAtomId() == "CA") or
									(a1.labelAtomId() == "CA" and a2.labelAtomId() == "O") or
									(a1.labelAtomId() == "N" and a2.labelAtomId() == "CB") or
									(a1.labelAtomId() == "CB" and a2.labelAtomId() == "N") or
									(a1.labelAtomId() == "C" and a2.labelAtomId() == "CB") or
									(a1.labelAtomId() == "CB" and a2.labelAtomId() == "C"))
								{
									minDist = 2.7;
								}
								break;

							case 2:
								if ((a1.labelAtomId() == "C" and a2.labelAtomId() == "N") or
									(a1.labelAtomId() == "N" and a2.labelAtomId() == "C"))
								{
									minDist = 2.7;
								}
								break;
						}
					}
				}
				catch (const exception& ex)
				{
					if (cif::VERBOSE)
						cerr << "err calculating nbc distance: " << ex.what() << endl;
					minDist = 2.8;
				}
			}
			
			mNonBondedContactRestraints.emplace_back(ra1, ra2, minDist, 0.02);
			nbc.insert(make_tuple(ra1, ra2));
			nbc.insert(make_tuple(ra2, ra1));
		}
	}

	// create reverse index (for dfcollector)
	mRef2AtomIndex = vector<size_t>(mReferencedAtoms.size(), kRefSentinel);
	for (size_t i = 0; i < mAtoms.size(); ++i)
	{
		AtomRef ar = ref(mAtoms[i]);
		assert(ar < mRef2AtomIndex.size());
		mRef2AtomIndex[ar] = i;
	}

	vector<pair<AtomRef,double>> densityAtoms;
	densityAtoms.reserve(mAtoms.size());
	
	transform(mAtoms.begin(), mAtoms.end(), back_inserter(densityAtoms),
		[this](const Atom& a) {
			double z = static_cast<int>(a.type());
			double weight = 1;
			double occupancy = a.occupancy();
			
			if (occupancy > 1)
				occupancy = 1;
			
			// TODO: cryo_em support
			
			return make_pair(ref(a), z * weight * occupancy);
		});

	mDensityRestraint.reset(new DensityRestraint(move(densityAtoms), xMap, mapWeight));
	
	// collect the restraints

	for (auto& r: mBondRestraints)				mRestraints.push_back(&r);
	for (auto& r: mAngleRestraints)				mRestraints.push_back(&r);
	for (auto& r: mTransPeptideRestraints)		mRestraints.push_back(&r);
	for (auto& r: mTorsionRestraints)			mRestraints.push_back(&r);
	for (auto& r: mPlanarityRestraints)			mRestraints.push_back(&r);
	for (auto& r: mChiralVolumeRestraints)		mRestraints.push_back(&r);
	for (auto& r: mNonBondedContactRestraints)	mRestraints.push_back(&r);
	if (mDensityRestraint)						mRestraints.push_back(mDensityRestraint.get());

	// report
	
	if (cif::VERBOSE)
		cout << "created " << mBondRestraints.size() << " bond restraints" << endl
			 << "created " << mAngleRestraints.size() << " angle restraints" << endl
			 << "created " << mTorsionRestraints.size() << " torsion restraints" << endl
			 << "created " << mPlanarityRestraints.size() << " plane restraints" << endl
			 << "created " << mTransPeptideRestraints.size() << " trans peptide restraints" << endl
			 << "created " << mChiralVolumeRestraints.size() << " chiral vol restraints" << endl
			 << "created " << mNonBondedContactRestraints.size() << " non-bonded-contact restraints" << endl
			 << endl;

	AtomLocationProvider loc(mReferencedAtoms);

	if (cif::VERBOSE > 1)
		for (auto r: mRestraints)
			r->print(loc);
}

AtomRef Minimizer::ref(const mmcif::Atom& atom)
{
	string atomID = atom.id();
	AtomRef result;
	
	auto k = mRefIndex.find(atomID);
	if (k != mRefIndex.end())
		result = k->second;
	else
	{
		result = mReferencedAtoms.size();
		mReferencedAtoms.push_back(atom);
		mRefIndex[atomID] = result;
	}
	
	return result;
}

void Minimizer::addLinkRestraints(const Monomer& a, const Monomer& b, const Link& link)
{
	auto& c1 = a.compound();
	auto& c2 = b.compound();
	
	auto getCompoundAtom = [&](const LinkAtom& la)
	{
		return la.compID == 1 ? c1.getAtomById(la.atomID) : c2.getAtomById(la.atomID);
	};
	
	auto getAtom = [&](const LinkAtom& la)
	{
		const Residue& r = la.compID == 1 ? a : b;
		return r.atomByID(la.atomID);
	};
	
	for (auto& b: link.bonds())
	{
		try
		{
			if (getCompoundAtom(b.atom[0]).typeSymbol == H or
				getCompoundAtom(b.atom[1]).typeSymbol == H)
			{
				continue;
			}
			
			Atom a1 = getAtom(b.atom[0]);
			Atom a2 = getAtom(b.atom[1]);

			mBondRestraints.emplace_back(ref(a1), ref(a2), b.distance, b.esd);
		}
		catch (const exception& ex)
		{
			if (cif::VERBOSE)
				cerr << "While processing bond restraints: " << ex.what() << endl;
			continue;
		}
	}

	for (auto& a: link.angles())
	{
		try
		{
			if (getCompoundAtom(a.atom[0]).typeSymbol == H or
				getCompoundAtom(a.atom[1]).typeSymbol == H or
				getCompoundAtom(a.atom[2]).typeSymbol == H)
			{
				continue;
			}
			
			Atom a1 = getAtom(a.atom[0]);
			Atom a2 = getAtom(a.atom[1]);
			Atom a3 = getAtom(a.atom[2]);

			mAngleRestraints.emplace_back(ref(a1), ref(a2), ref(a3), a.angle, a.esd);
		}
		catch (const exception& ex)
		{
			if (cif::VERBOSE)
				cerr << "While processing angle restraints: " << ex.what() << endl;
			continue;
		}
	}
	
//	for (auto& a: link.torsions())
//	{
//		if (a.esd == 0)
//			continue;
//		
//		try
//		{
//			if (getCompoundAtom(a.atom[0]).typeSymbol == H or
//				getCompoundAtom(a.atom[1]).typeSymbol == H or
//				getCompoundAtom(a.atom[2]).typeSymbol == H or
//				getCompoundAtom(a.atom[3]).typeSymbol == H)
//			{
//				continue;
//			}
//			
//			Atom a1 = getAtom(a.atom[0]);
//			Atom a2 = getAtom(a.atom[1]);
//			Atom a3 = getAtom(a.atom[2]);
//			Atom a4 = getAtom(a.atom[3]);
//
//			mTorsionRestraints.emplace_back(ref(a1), ref(a2), ref(a3), ref(a4), a.angle, a.esd, a.period);
//		}
//		catch (const exception& ex)
//		{
//			if (cif::VERBOSE)
//				cerr << "While processing torsion restraints: " << ex.what() << endl;
//			continue;
//		}
//	}
	
	for (auto& cv: link.chiralCentres())
	{
		try
		{
			if (getCompoundAtom(cv.atom[0]).typeSymbol == H or
				getCompoundAtom(cv.atom[1]).typeSymbol == H or
				getCompoundAtom(cv.atom[2]).typeSymbol == H)
			{
				continue;
			}
			
			Atom cc = getAtom(cv.atomCentre);
			Atom a1 = getAtom(cv.atom[0]);
			Atom a2 = getAtom(cv.atom[1]);
			Atom a3 = getAtom(cv.atom[2]);

			auto volume = link.chiralVolume(cv.id);

			mChiralVolumeRestraints.emplace_back(ref(cc), ref(a1), ref(a2), ref(a3), volume);
		}
		catch (const exception& ex)
		{
			if (cif::VERBOSE)
				cerr << "While processing chiral volume restraints: " << ex.what() << endl;
			continue;
		}
	}
	
	for (auto& p: link.planes())
	{
		try
		{
			vector<AtomRef> atoms;
			
			for (auto a: p.atoms)
			{
				if (getCompoundAtom(a).typeSymbol == H)
					continue;
				
				atoms.push_back(ref(getAtom(a)));
			}
			
			if (atoms.size() > 3)
				mPlanarityRestraints.emplace_back(PlanarityRestraint{move(atoms), p.esd});
		}
		catch (const exception& ex)
		{
			if (cif::VERBOSE)
				cerr << "While processing planarity restraints: " << ex.what() << endl;
			continue;
		}
	}
}

void Minimizer::printStats()
{
	AtomLocationProvider loc(mReferencedAtoms);

	double bondScore = rmsz(loc, mBondRestraints);
	double angleScore = rmsz(loc, mAngleRestraints);
	double torsionScore = rmsz(loc, mTorsionRestraints);
	double chiralityVolumeScore = rmsz(loc, mChiralVolumeRestraints);
	double planarityScore = rmsz(loc, mPlanarityRestraints);
	double transpeptideScore = rmsz(loc, mTransPeptideRestraints);
	double nbcScore = rmsz(loc, mNonBondedContactRestraints);
	double densityScore = mDensityRestraint->f(loc);

	cerr << "  Bonds:              " << bondScore << endl
		 << "  Angles:             " << angleScore << endl
		 << "  Torsion:            " << torsionScore << endl
		 << "  Chirality:          " << chiralityVolumeScore << endl
		 << "  Planarity:          " << planarityScore << endl
		 << "  Transpeptide:       " << transpeptideScore << endl
		 << "  Non-Bonded-Contact: " << nbcScore << endl
		 << "  Density:            " << densityScore << endl;
}

double Minimizer::score()
{
	AtomLocationProvider loc(mReferencedAtoms);
	return score(loc);
}

double Minimizer::score(const AtomLocationProvider& loc)
{
	double result = 0;
	for (auto r: mRestraints)
		result += r->f(loc);
	
	if (cif::VERBOSE > 3)
		cout << "score: " << result << endl;

	return result;
}

// --------------------------------------------------------------------

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h> // for debugging norm of gradient

class GSLAtomLocation : public AtomLocationProvider
{
  public:
	GSLAtomLocation(vector<Atom>& atoms, const vector<Point>& fixedAtoms,
		const vector<size_t>& index, const gsl_vector* v)
		: AtomLocationProvider(atoms), mFixedLocations(fixedAtoms), mIndex(index), mV(v)
	{
		assert(mIndex.size() == mFixedLocations.size());

		if (cif::VERBOSE > 2)
		{
			for (size_t i = 0; i < mIndex.size(); ++i)
			{
				size_t ri = mIndex[i];
				if (ri == kRefSentinel)
					continue;
		
				DPoint p = {
					gsl_vector_get(mV, ri * 3    ),
					gsl_vector_get(mV, ri * 3 + 1),
					gsl_vector_get(mV, ri * 3 + 2)
				};
		
				cout << mAtoms[i] << p << endl;
			}
		}
	}
	
	virtual DPoint operator[](AtomRef atom) const;

	void storeLocations();

  private:
	const vector<Point>& mFixedLocations;
	const vector<size_t>& mIndex;
	const gsl_vector* mV;
};

DPoint GSLAtomLocation::operator[](AtomRef atomID) const
{
	assert(atomID < mIndex.size());
	
	size_t ix = mIndex.at(atomID);

	if (ix == kRefSentinel)
		return mFixedLocations.at(atomID);
	
	return DPoint(
		gsl_vector_get(mV, ix * 3 + 0),
		gsl_vector_get(mV, ix * 3 + 1),
		gsl_vector_get(mV, ix * 3 + 2)
	);
}

void GSLAtomLocation::storeLocations()
{
	for (size_t i = 0; i < mIndex.size(); ++i)
	{
		size_t ri = mIndex[i];
		if (ri == kRefSentinel)
			continue;

		DPoint p = {
			gsl_vector_get(mV, ri * 3    ),
			gsl_vector_get(mV, ri * 3 + 1),
			gsl_vector_get(mV, ri * 3 + 2)
		};

		mAtoms[i].location(p);
	}
}

// --------------------------------------------------------------------

class GSLDFCollector : public DFCollector
{
  public:
	GSLDFCollector(const vector<mmcif::Atom>& atoms, const vector<size_t>& index, gsl_vector* df)
		: mAtoms(atoms), mIndex(index), mDF(df)
	{
		for (size_t ix: mIndex)
		{
			if (ix == kRefSentinel)
				continue;
			
			gsl_vector_set(mDF, ix * 3 + 0, 0.0);
			gsl_vector_set(mDF, ix * 3 + 1, 0.0);
			gsl_vector_set(mDF, ix * 3 + 2, 0.0);
		}
	}

	~GSLDFCollector();
	
	virtual void add(AtomRef atom, double dx, double dy, double dz);

  private:

	// for debugging
	string label(AtomRef atom) const
	{
		string atomName = " " + mAtoms[atom].labelAtomId();
		atomName += string(5 - atomName.length(), ' ');
		return to_string(mAtoms[atom].labelSeqId()) + atomName;
	}

	const vector<mmcif::Atom>& mAtoms;
	const vector<size_t>& mIndex;
	gsl_vector* mDF;
};

GSLDFCollector::~GSLDFCollector()
{
	if (cif::VERBOSE > 1)
	{
		cerr << string(cif::get_terminal_width(), '-') << endl
			 << "Collected gradient: " << endl;
		
		for (size_t i = 0; i < mAtoms.size(); ++i)
		{
			size_t ix = mIndex[i];
			if (ix == kRefSentinel)
				continue;
			
			double dx = gsl_vector_get(mDF, ix * 3 + 0);
			double dy = gsl_vector_get(mDF, ix * 3 + 1);
			double dz = gsl_vector_get(mDF, ix * 3 + 2);
			
			cerr << "atom: " << label(i) << " d: " << setprecision(10) << dx << " " << dy << " " << dz << endl;
		}

		cerr << string(cif::get_terminal_width(), '-') << endl;
	}
}

void GSLDFCollector::add(AtomRef atom, double dx, double dy, double dz)
{
	assert(atom < mIndex.size());

	size_t ix = mIndex[atom];
	if (ix != kRefSentinel)
	{
		gsl_vector_set(mDF, ix * 3 + 0, gsl_vector_get(mDF, ix * 3 + 0) + dx);
		gsl_vector_set(mDF, ix * 3 + 1, gsl_vector_get(mDF, ix * 3 + 1) + dy);
		gsl_vector_set(mDF, ix * 3 + 2, gsl_vector_get(mDF, ix * 3 + 2) + dz);
		
		if (cif::VERBOSE > 1)
			cerr << "atom: " << label(atom) << " d: " << setprecision(10) << dx << ", " << dy << ", " << dz << endl;
	}
}

// --------------------------------------------------------------------

class GSLMinimizer : public Minimizer
{
  public:
	GSLMinimizer(const Polymer& poly, int first, int last, const BondMap& bonds,
		const XMap& xMap, float mapWeight, float plane5AtomsESD)
		: Minimizer(poly, first, last, bonds, xMap, mapWeight, plane5AtomsESD)
	{
		for (auto& a: mReferencedAtoms)
			mFixedLocations.push_back(a.location());
	}
	
	~GSLMinimizer()
	{
		if (m_s != nullptr)
			gsl_multimin_fdfminimizer_free(m_s);
	}
	
	virtual double refine(bool storeAtoms);
	virtual vector<pair<string,mmcif::Point>> getAtoms() const;
	virtual void storeAtomLocations();

  private:

	static double F(const gsl_vector* v, void* params);
	static void Df(const gsl_vector* v, void* params, gsl_vector* df);
	static void Fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df);

	double F(const gsl_vector* v);
	void Df(const gsl_vector* v, gsl_vector* df);
	void Fdf (const gsl_vector *x, double *f, gsl_vector *df);

	vector<Point> mFixedLocations;
	gsl_multimin_fdfminimizer* m_s = nullptr;
};

double GSLMinimizer::refine(bool storeAtoms)
{
	const size_t iterations = 4000;
	
	gsl_multimin_function_fdf fdf = {};
	fdf.f = &GSLMinimizer::F;
	fdf.df = &GSLMinimizer::Df;
	fdf.fdf = &GSLMinimizer::Fdf;
	fdf.n = mAtoms.size() * 3;
	fdf.params = this;

	auto T = gsl_multimin_fdfminimizer_conjugate_pr;
	auto x = gsl_vector_alloc(3 * mAtoms.size());

	size_t ix = 0;
	for (auto& a: mAtoms)
	{
		auto l = a.location();
		gsl_vector_set(x, ix++, l.mX);
		gsl_vector_set(x, ix++, l.mY);
		gsl_vector_set(x, ix++, l.mZ);
	}
	
	m_s = gsl_multimin_fdfminimizer_alloc(T, 3 * mAtoms.size());

	float tolerance = 0.06f;
	double stepSize = 0.1 * gsl_blas_dnrm2(x);

	gsl_multimin_fdfminimizer_set(m_s, &fdf, x, stepSize, tolerance);
	
	double gradLim = sqrt(mRestraints.size()) * 0.15;
	if (gradLim < 0.3)
		gradLim = 0.3;

	for (size_t i = 0; i < iterations; ++i)
	{
		int status = gsl_multimin_fdfminimizer_iterate(m_s);
		
		if (status != 0)
		{
			if (status != GSL_ENOPROG)
				cerr << "Unexpected result from gsl_multimin_fdfminimizer_iterate: " << status << endl;
			else if (cif::VERBOSE)
				cerr << "Minimizer stopped at iteration " << i << " at " << m_s->f << endl;
			break;
		}

		status = gsl_multimin_test_gradient(m_s->gradient, gradLim);

		if (cif::VERBOSE > 1)
		{
			double norm = gsl_blas_dnrm2(m_s->gradient);
			cout << "iteration number " << i << " with f: " << m_s->f
			      << " status from gsl_multimin_test_gradient() " << status << " for norm "
			      << norm << endl;
		}
		
		if (status == GSL_SUCCESS)
		{
			if (cif::VERBOSE)
				cerr << "Minimum found at iteration " << i << " at " << m_s->f << endl;
			break;
		}
		
		if (status != GSL_CONTINUE)
			break;
	}
	
	gsl_vector_free(x);
	
	if (storeAtoms)
		storeAtomLocations();
	
	return m_s->f;
}

vector<pair<string,mmcif::Point>> GSLMinimizer::getAtoms() const
{
	vector<pair<string,mmcif::Point>> result;

	for (size_t i = 0; i < mRef2AtomIndex.size(); ++i)
	{
		size_t ri = mRef2AtomIndex[i];
		if (ri == kRefSentinel)
			continue;

		DPoint p = {
			gsl_vector_get(m_s->x, ri * 3    ),
			gsl_vector_get(m_s->x, ri * 3 + 1),
			gsl_vector_get(m_s->x, ri * 3 + 2)
		};
		
		result.emplace_back(mReferencedAtoms[i].id(), p);
	}
	
	return result;
}

void GSLMinimizer::storeAtomLocations()
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, m_s->x);
	loc.storeLocations();
}

double GSLMinimizer::F(const gsl_vector* v, void* params)
{
	GSLMinimizer* self = reinterpret_cast<GSLMinimizer*>(params);
	return self->F(v);
}

void GSLMinimizer::Df(const gsl_vector* v, void* params, gsl_vector* df)
{
	GSLMinimizer* self = reinterpret_cast<GSLMinimizer*>(params);
	self->Df(v, df);
}

void GSLMinimizer::Fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df)
{
	GSLMinimizer* self = reinterpret_cast<GSLMinimizer*>(params);
	self->Fdf(v, f, df);

	if (cif::VERBOSE > 1)
		cout << "FDF => " << setprecision(10) << *f << endl;
}

double GSLMinimizer::F(const gsl_vector* v)
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, v);
	
//	return score(loc);
	auto F = score(loc);
	if (cif::VERBOSE > 1)
		cout << "F => " << setprecision(10) << F << endl;
	return F;
}

void GSLMinimizer::Df(const gsl_vector* v, gsl_vector* df)
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, v);
	
	GSLDFCollector c(mReferencedAtoms, mRef2AtomIndex, df);
	
	for (auto r: mRestraints)
		r->df(loc, c);
}

void GSLMinimizer::Fdf (const gsl_vector *x, double *f, gsl_vector *df)
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, x);
	*f = score(loc);

	GSLDFCollector c(mReferencedAtoms, mRef2AtomIndex, df);
	for (auto r: mRestraints)
		r->df(loc, c);
}

// --------------------------------------------------------------------

Minimizer* Minimizer::create(const string& algorithm,	
	const Polymer& poly, int first, int last, const BondMap& bonds,
	const XMap& xMap, float mapWeight, float plane5AtomsESD)
{
	Minimizer* result = nullptr;

	if (algorithm == "gsl")
		result = new GSLMinimizer(poly, first, last, bonds, xMap, mapWeight, plane5AtomsESD);
	else
		throw runtime_error("Unknown algorithm: " + algorithm);
	
	return result;
}
