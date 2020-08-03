/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 19 februari, 2018
*/

// test 3fvl

#include "pdb-redo.h"

#include <iomanip>
#include <numeric>
#include <future>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>

#include "cif++/Secondary.hpp"
#include "cif++/CifUtils.hpp"

#include "minimizer.h"
#include "ramachandran.h"

using namespace std;

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace ba = boost::algorithm;

using mmcif::Atom;
using mmcif::Point;
using mmcif::Structure;
using mmcif::Monomer;
using mmcif::BondMap;

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Coord_map;
using clipper::Coord_frac;

// --------------------------------------------------------------------

void Compare(Structure& a, Structure& b)
{
	vector<Atom> atomsA = a.atoms(), atomsB = b.atoms();

	if (atomsA.size() != atomsB.size())
		throw runtime_error("Ongelijk aantal atomen");

	vector<double> dists;
	
	for (size_t i = 0; i < atomsA.size(); ++i)
	{
		Atom aa = atomsA[i];
		Atom ab = atomsB[i];
		
		double d = Distance(aa, ab);
		if (d == 0)
			continue;
		
		dists.push_back(d);
	}
	
	double median = 0;
	if (dists.size() % 2 == 1)
		median = dists[dists.size() / 2];
	else if (not dists.empty())
		median = (dists[dists.size() / 2 - 1] + dists[dists.size() / 2]) / 2;
	
	cout << "Number of differing atom positions: " << dists.size() << endl
		 << "Average distance: " << accumulate(dists.begin(), dists.end(), 0.0) / dists.size() << endl
		 << "Median distance: " << median << endl;
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("xyzin-1",				po::value<string>(),	"coordinates file")
		("xyzin-2",				po::value<string>(),	"coordinates file")

		("help,h",										"Display help message")
		("version",										"Print version")

		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)")
		;

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin-1", 1);
	p.add("xyzin-2", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	
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
	
	if (vm.count("xyzin-1") == 0 or vm.count("xyzin-2") == 0)
	{
		cerr << "Input files not specified" << endl;
		exit(1);
	}

	mmcif::File f1(vm["xyzin-1"].as<string>());
	Structure s1(f1);
	
	mmcif::File f2(vm["xyzin-2"].as<string>());
	Structure s2(f2);

	Compare(s1, s2);
	
	return 0;
}
	
