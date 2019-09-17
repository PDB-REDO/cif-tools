/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 19 februari, 2018
*/

#include "pdb-redo.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <ctgmath>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>

#include <zeep/xml/document.hpp>
#include <zeep/xml/serialize.hpp>
#include <zeep/xml/writer.hpp>

#include "cif++/Cif++.h"
#include "cif++/Structure.h"
#include "cif++/Compound.h"
#include "cif++/AtomShape.h"
#include "cif++/CifUtils.h"
#include "cif++/BondMap.h"
#include "cif++/MapMaker.h"
#include "cif++/ResolutionCalculator.h"

#include "HBondTraits.h"

#include "cif++/mrsrc.h"
#include "svm++.h"

using namespace std;
using namespace clipper;
using data32::Flag;
using data32::F_phi;
using data32::F_sigF;
using data32::Phi_fom;
using mmcif::kPI;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace zx = zeep::xml;
namespace c = mmcif;

// --------------------------------------------------------------------

class ChironHelper
{
  public:
	ChironHelper(const mmcif::Compound& compound, const string& centre)
		: mCompound(compound), mCentre(centre)
		, mAtoms(mCompound.atoms()), mBonds(mCompound.bonds())
	{
	}
	
	bool canSwapChains(vector<tuple<string,string>>& swapAtoms);
	bool canSwapIsomers(const mmcif::Residue& res, string& swapCompound, vector<tuple<string,string>>& swapAtoms);
	bool canPullThroughPlane(const mmcif::Residue& r, mmcif::Point& newCoordinates);
	
  private:

	bool canSwapChains(const string& atom1, const string& atom2,
		vector<tuple<string,string>>& swapAtoms);

	bool canSwapSubChain(string atom1, string atom2,
		set<string> visited1, set<string> visited2,
		vector<tuple<string,string>>& swapAtoms);

	const mmcif::CompoundAtom& atomForID(const string& id) const
	{
		auto i = find_if(mAtoms.begin(), mAtoms.end(),
			[&](const mmcif::CompoundAtom& a) { return a.id == id; });
		if (i == mAtoms.end())
			throw runtime_error("Could not find atom " + id + " in compound " + mCompound.id());
		return *i;
	}

	size_t countChiralErrors(const mmcif::Residue& res, const mmcif::Compound& c) const
	{
		vector<tuple<string,string>> remapped;
		return countChiralErrors(res, c, remapped);
	}
	
	size_t countChiralErrors(const mmcif::Residue& res, const mmcif::Compound& c, const vector<tuple<string,string>>& remapped) const;

	const mmcif::Compound&			mCompound;
	string							mCentre;
	vector<mmcif::CompoundAtom>	mAtoms;
	vector<mmcif::CompoundBond>	mBonds;
};

bool ChironHelper::canPullThroughPlane(const mmcif::Residue& r, mmcif::Point& newCoordinates)
{
	bool result = false;

	// see if one of the bonds is to a hydrogen
	vector<string> l;
	for (auto& b: mBonds)
	{
		string linked;
		
		if (b.atomID[0] == mCentre)
			linked = b.atomID[1];
		else if (b.atomID[1] == mCentre)
			linked = b.atomID[0];
		
		if (not linked.empty() and atomForID(linked).typeSymbol != mmcif::H)
			l.push_back(linked);
	}
	
	// Only when we have three atoms other than H bonded we can pull
	// the centre through the plane
	if (l.size() == 3)
	{
		result = true;
		
		// Construct plane based on three points
		// a.x + b.y + c.z + d = 0

		auto centre = r.atomByID(mCentre).location();
		
		// Normal vector on the plane is:
		auto p1 = r.atomByID(l[0]).location();
		auto p2 = r.atomByID(l[1]).location();
		auto p3 = r.atomByID(l[2]).location();
		
		auto n = CrossProduct(p2 - p1, p3 - p1);
		
		float a = n.getX();
		float b = n.getY();
		float c = n.getZ();
		float d = p1.getX() * a + p1.getY() * b + p1.getZ() * c;
		
		d -= a * centre.getX() + b * centre.getY() + c * centre.getZ();
		
		// point closest to centre for this plane:
		mmcif::Point o {
			a * d / (a * a + b * b + c * c),
			b * d / (a * a + b * b + c * c),
			c * d / (a * a + b * b + c * c)
		};
		o += centre;
		
		auto l = o - centre;
		
		newCoordinates = o + l / 100.f;
	}
	
	return result;
}

bool ChironHelper::canSwapChains(vector<tuple<string,string>>& swapAtoms)
{
	vector<string> l;
	
	for (auto& b: mBonds)
	{
		if (b.atomID[0] == mCentre)
			l.push_back(b.atomID[1]);
		else if (b.atomID[1] == mCentre)
			l.push_back(b.atomID[0]);
	}
	
	bool result = false;

	for (size_t i = 0; i + 1 < l.size(); ++i)
	{
		for (size_t j = i + 1; j < l.size(); ++j)
		{
			if (canSwapChains(l[i], l[j], swapAtoms))
			{
				result = true;
				break;
			}
		}
	}

	return result;
}

bool ChironHelper::canSwapChains(const string& atom1, const string& atom2,
		vector<tuple<string,string>>& swapAtoms)
{
	set<string> v1, v2;
	
	v1.insert(mCentre);
	v2.insert(mCentre);
	
	return canSwapSubChain(atom1, atom2, v1, v2, swapAtoms);
}

bool ChironHelper::canSwapIsomers(const mmcif::Residue& res, string& swapCompound, vector<tuple<string,string>>& swapAtoms)
{
	bool result = false;
	
	if (VERBOSE or not mCompound.isSugar())
	{
		auto isomers = mCompound.isomers();
	
		if (not isomers.empty())
		{
			auto chiralErrInitial = countChiralErrors(res, mCompound);
			auto chiralErr = chiralErrInitial;
			
			if (VERBOSE > 1)
				cerr << "Trying to swap isomers, initial error count is " << chiralErrInitial << endl;
			
			for (auto i: isomers)
			{
				auto c = mmcif::Compound::create(i);
				
				vector<tuple<string,string>> m = c->mapToIsomer(mCompound);
				
				if (VERBOSE > 2)
				{
					for (auto& a: m)
						cerr << "  " << get<0>(a) << " => " << get<1>(a) << endl;
				}
				
				auto err = countChiralErrors(res, *c, m);
	
				if (VERBOSE > 1)
					cerr << "isomer " << i << " has " << err << " errors" << endl;
	
				if (chiralErr > err)
				{
					chiralErr = err;
					result = true;
					swapCompound = i;
					swap(swapAtoms, m);
					
					if (VERBOSE > 1)
						cerr << "Err count decreased to " << chiralErr << endl; 
	
					if (chiralErr == 0)	 // we're done.
						break;
				}
			}
	
			if (result and mCompound.isSugar())
			{
				if (VERBOSE)
					cerr << "Since residue " << res.compoundID() << " is a sugar, it will not be swapped with " << swapCompound << endl;
				result = false;
			}
		}
	}
	
	return result;
}

size_t ChironHelper::countChiralErrors(const mmcif::Residue& res, const mmcif::Compound& c, const vector<tuple<string,string>>& mapping) const
{
	size_t result = 0;
	
	for (auto cc: c.chiralCentres())
	{
		if (cc.volumeSign == mmcif::both)
			continue;
		
		auto rename = [&](const string& name) -> string
		{
			string result;
			
			auto i = find_if(mapping.begin(), mapping.end(), [&](auto& m) { return get<0>(m) == name; });
			if (i == mapping.end())
			{
//				if (VERBOSE > 1 and not mapping.empty())
//					cerr << "no mapping found for atom " << name << " in " << c.id() << endl;
				 result = name;
			}
			else
				result = get<1>(*i);

			return result;
		};
		
		try
		{
			auto centre = res.atomByID(rename(cc.atomIDCentre));
			auto atom1 = res.atomByID(rename(cc.atomID[0]));
			auto atom2 = res.atomByID(rename(cc.atomID[1]));
			auto atom3 = res.atomByID(rename(cc.atomID[2]));

			auto chiralVolume = DotProduct(atom1.location() - centre.location(),
				CrossProduct(atom2.location() - centre.location(), atom3.location() - centre.location()));

			if ((chiralVolume < 0 and cc.volumeSign == mmcif::positiv) or
				(chiralVolume > 0 and cc.volumeSign == mmcif::negativ))
			{
				if (VERBOSE > 1)
					cerr << "chiral error in " << c.id() << " around " << cc.atomIDCentre << " with volume: " << chiralVolume << endl;
				
				++result;
			}
		}
		catch (const exception& ex)
		{
			if (VERBOSE)
				cerr << "Missing atom in counting chiral errors: " << ex.what() << endl;
		}
	}
	
	return result;
}

bool ChironHelper::canSwapSubChain(string atom1, string atom2,
		set<string> visited1, set<string> visited2,
		vector<tuple<string,string>>& swapAtoms)
{
	auto& a1 = atomForID(atom1);
	auto& a2 = atomForID(atom2);
	
	bool result = false;
	for (;;)
	{
		if (a1.typeSymbol != a2.typeSymbol)
			break;
		
		if (a1.typeSymbol == mmcif::H)
		{
			result = true;
			break;
		}

		vector<pair<string,mmcif::BondType>> l1, l2;
		
		visited1.insert(atom1);
		visited2.insert(atom2);
		
		for (auto& b: mBonds)
		{
			if (atom1 == b.atomID[0] and not visited1.count(b.atomID[1]))
				l1.push_back(make_pair(b.atomID[1], b.type));
			else if (atom1 == b.atomID[1] and not visited1.count(b.atomID[0]))
				l1.push_back(make_pair(b.atomID[0], b.type));

			if (atom2 == b.atomID[0] and not visited2.count(b.atomID[1]))
				l2.push_back(make_pair(b.atomID[1], b.type));
			else if (atom2 == b.atomID[1] and not visited2.count(b.atomID[0]))
				l2.push_back(make_pair(b.atomID[0], b.type));
		}
		
		if (l1.size() != l2.size())
			break;
		
		vector<tuple<string,string>> subSwap;		
		
		auto test = [&](int a, int b) -> bool
		{
			bool result = l1[a].second == l2[b].second and
				(l1[a].first == l2[a].first or canSwapSubChain(l1[a].first, l2[b].first, visited1, visited2, subSwap));
			if (not result)
				subSwap.clear();
			return result;
		};
		
		auto add = [&](initializer_list<pair<int,int>> l)
		{
			if (atom1 != atom2)
				swapAtoms.emplace_back(atom1, atom2);
			swapAtoms.insert(swapAtoms.end(), subSwap.begin(), subSwap.end());
			result = true;
		};
		
		switch (l1.size())
		{
			case 0:
				add({});
				break;
			
			case 1:
				if (test(0, 0))
					add({ make_pair(0, 0) });
				break;
			
			case 2:
				if (test(0, 0) and test(1, 1))
					add({ make_pair(0, 0), make_pair(1, 1) });
				else if (test(0, 1) and test(1, 0))
					add({ make_pair(0, 1), make_pair(1, 0) });
				break;
			
			case 3:
				if (test(0, 0))
				{
					if (test(1, 1) and test(2, 2)) 			add({ make_pair(0, 0), make_pair(1, 1), make_pair(2, 2) });
					else if (test(1, 2) and test(2, 1))		add({ make_pair(0, 0), make_pair(1, 2), make_pair(2, 1) });
				}
				else if (test(0, 1))
				{
					if (test(1, 2) and test(2, 0))			add({ make_pair(0, 1), make_pair(1, 2), make_pair(2, 0) });
					else if (test(1, 0) and test(2, 2))		add({ make_pair(0, 1), make_pair(1, 0), make_pair(2, 2) });
				}
				else if (test(0, 2))
				{
					if (test(1, 0) and test(2, 1))			add({ make_pair(0, 2), make_pair(1, 0), make_pair(2, 1) });
					else if (test(1, 1) and test(2, 0))		add({ make_pair(0, 2), make_pair(1, 1), make_pair(2, 0) });
				}
	
				break;
			
			default:
				throw runtime_error("unimplemented number of bonds");
		}

		break;
	}

	if (VERBOSE > 2)
		cerr << "canSwap(" << atom1 << ", " << atom2 << ") => " << boolalpha << result << endl;

	return result;
}

// --------------------------------------------------------------------

int Process(mmcif::Structure& structure, const mmcif::Residue& res, const mmcif::Compound& compound)
{
	int result = 0;
	
	if (VERBOSE > 1)
		cerr << "Process " << res.compoundID() << " " << res.asymID() << res.seqID() << endl;
	
	for (auto cc: compound.chiralCentres())
	{
		if (cc.volumeSign == mmcif::both)
			continue;
		
		try
		{
			auto centre = res.atomByID(cc.atomIDCentre);
			auto atom1 = res.atomByID(cc.atomID[0]);
			auto atom2 = res.atomByID(cc.atomID[1]);
			auto atom3 = res.atomByID(cc.atomID[2]);

			auto chiralVolume = DotProduct(atom1.location() - centre.location(),
				CrossProduct(atom2.location() - centre.location(), atom3.location() - centre.location()));

			if (VERBOSE)
			{
				cerr << "chiral volume for " << res.compoundID() << " " << res.asymID() << res.seqID()
					 << " with centre " << cc.atomIDCentre
					 << " is " << chiralVolume << " and should be "
					 << (cc.volumeSign == mmcif::positiv ? "positive" : "negative") << endl;
			}

			if ((chiralVolume < 0 and cc.volumeSign == mmcif::positiv) or
				(chiralVolume > 0 and cc.volumeSign == mmcif::negativ))
			{
				cerr << "Error in chiral volume for " << res.compoundID() << " " << res.asymID() << res.seqID()
					 << " with centre " << cc.atomIDCentre;
				
				if (chiralVolume < 0)
					cerr << " volume should be positive but is negative: " << chiralVolume << endl;
				else
					cerr << " volume should be negative but is positive: " << chiralVolume << endl;

				ChironHelper test(compound, cc.atomIDCentre);
				
				vector<tuple<string,string>> swapAtoms;
				string swapCompound;
				mmcif::Point newCoordinates;
				
				if (test.canSwapChains(swapAtoms))
				{
					++result;
					
					cerr << "Flipping labels: ";
					
					for (auto p: swapAtoms)
					{
						string na1, na2;
						tie(na1, na2) = p;
						
						cerr << "{ " << na1 << " and " << na2 << " }, ";

						auto a1 = res.atomByID(na1);
						auto a2 = res.atomByID(na2);

						structure.swapAtoms(a1, a2);
					}

					cerr << endl;
				}
				else if (test.canSwapIsomers(res, swapCompound, swapAtoms))
				{
					++result;
					cerr << "Replacing with isomer " << swapCompound;

					if (swapAtoms.empty())
						cerr << endl;
					else
					{
						cerr << ", swapping atom labels: ";

						for (auto p: swapAtoms)
						{
							string na1, na2;
							tie(na1, na2) = p;
							
							cerr << "{ " << na1 << " -> " << na2 << " }, ";
						}

						cerr << endl;
					}
					
					structure.changeResidue(res, swapCompound, swapAtoms);
				}
				else if (test.canPullThroughPlane(res, newCoordinates))
				{
					++result;

					structure.moveAtom(centre, newCoordinates);
					
					chiralVolume = DotProduct(atom1.location() - newCoordinates,
						CrossProduct(atom2.location() - newCoordinates, atom3.location() - newCoordinates));
		
					cerr << "Pulling it through the plane from " << centre.location() << " to " << newCoordinates << endl
						 << "  new chiral volume is " << chiralVolume << endl;
				}
				else
					cerr << "Cannot fix this problem" << endl;
			}
		}
		catch (const runtime_error& ex)
		{
			if (VERBOSE)
				cerr << ex.what() << endl;
			
			continue;
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

int Process(mmcif::Structure& structure)
{
	int result = 0;
	
	for (auto& res: structure.nonPolymers())
	{
		// skip waters...
		if (res.isWater())
			continue;
		
		auto& compound = res.compound();

		if (compound.chiralCentres().empty())
			continue;
		
		result += Process(structure, res, compound);
	}

	for (auto& poly: structure.polymers())
	{
		for (auto& m: poly)
		{
			auto& compound = m.compound();

			if (compound.chiralCentres().empty())
				continue;
			
			result += Process(structure, m, compound);
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("xyzin",				po::value<string>(),	"coordinates file")
		("output,o",			po::value<string>(),	"Write output to this file instead of stdout")
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("xyzin") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}
	
	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();

	mmcif::File f(vm["xyzin"].as<string>());
	mmcif::Structure structure(f);

	int result = Process(structure);

	if (result)
	{
		if (vm.count("output"))
			structure.getFile().save(vm["output"].as<string>());
		else
			structure.getFile().file().save(cout);
	}
	
	return result;
}
