/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 19 februari, 2018
*/

// test 3fvl

#include "pdb-redo.h"

#include <iomanip>
#include <fstream>
#include <filesystem>

#include <boost/program_options.hpp>

#include "cif++/Statistics.hpp"
#include "cif++/CifUtils.hpp"

#include "minimizer.h"

using namespace std;

namespace po = boost::program_options;
namespace fs = std::filesystem;

using mmcif::Atom;
using mmcif::Point;
using mmcif::Structure;
using mmcif::Monomer;
using mmcif::BondMap;

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Coord_map;
using clipper::Coord_frac;

typedef mmcif::MapMaker<float> MapMaker;

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("hklin",				po::value<string>(),	"reflections file")
		("xyzin",				po::value<string>(),	"coordinates file")
		("xyzout,o",			po::value<string>(),	"output coordinates to this file")
		("asym-id",				po::value<string>(),	"Asymetric unit to refine")
		("res-first",			po::value<int>(),		"Sequence number for first residue to refine, default = 1")
		("res-last",			po::value<int>(),		"Sequence number for last residue to refine, default is last in sequence")
		("iterations,i",		po::value<uint32_t>(),	"Maximum number of iterations, default is nr of moving atoms x 1000")
		("algorithm",			po::value<string>(),	"Optimisation algorithm (either gsl or sa)")
		
		("weight-density",		po::value<float>(),		"Weight for density score")
//		("nudge-offset",		po::value<float>(),		"Maximum offset for nudging atom in simulated annealing") 
		("sampling-rate",		po::value<float>(),		"Sampling rate")
		("plane-5-atoms-esd",	po::value<float>(),		"ESD for the atoms in the 5 atom peptide bond plane, default is 0.11")

		("stats",										"Calculated statistics for refined residues")

		("dict",				po::value<vector<string>>(),
														"Dictionary file containing restraints for residues in this specific target, can be specified multiple times.")
		
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
	p.add("hklin", 1);
	p.add("xyzin", 1);
	p.add("xyzout", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	
	if (fs::exists("re-ref.cfg"))
	{
		std::ifstream cfgFile("re-ref.cfg");
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}
	
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
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

	if (vm.count("asym-id") == 0)
	{
		cerr << "Asym-id not specified" << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	if (vm.count("dict"))
	{
		for (auto dict: vm["dict"].as<vector<string>>())
			mmcif::CompoundFactory::instance().pushDictionary(dict);
	}

	string algo = "gsl";
	if (vm.count("algorithm"))
		algo = vm["algorithm"].as<string>();

	float mapWeight = 60;
	if (vm.count("weight-density"))
		mapWeight = vm["weight-density"].as<float>();

	mmcif::File f(vm["xyzin"].as<string>());
	Structure structure(f);
	
	fs::path mtzFile = vm["hklin"].as<string>();

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

	float plane5AtomsESD = 0.11;
	if (vm.count("plane-5-atoms-esd"))
		plane5AtomsESD = vm["plane-5-atoms-esd"].as<float>();

	mmcif::MapMaker<float> mm;
	mm.loadMTZ(mtzFile, samplingRate);

	mmcif::BondMap bm(structure);

	unique_ptr<Minimizer> minimizer;

	string asymID = vm["asym-id"].as<string>();
	
	int resFirst = 1;
	if (vm.count("res-first"))
		resFirst = vm["res-first"].as<int>();
	
	int resLast = numeric_limits<int>::max();
	if (vm.count("res-last"))
		resLast = vm["res-last"].as<int>();

	if (vm.count("stats"))
	{
		mmcif::StatsCollector collector(mm, structure, false /*electronScattering*/);
		auto r = collector.collect(asymID, resFirst, resLast);
	
		cout << "RESIDUE" << '\t'
			 << "RSR" << '\t'
			 << "SRSR" << '\t'
			 << "RSCCS" << '\t'
			 << "NGRID" << '\t'
			 << endl;
	
		for (auto i: r)
		{
			cout << fixed << setprecision(3)
				 << i.compID << '_' << i.asymID << '_' << i.seqID << '\t'
				 << i.RSR << '\t'
				 << i.SRSR << '\t'
				 << i.RSCCS << '\t'
				 << i.ngrid << '\t'
				 << endl;
		}
	}
	
	for (auto& poly: structure.polymers())
	{
		if (poly.asymID() != asymID)
			continue;
			
		minimizer.reset(Minimizer::create(algo, poly, resFirst, resLast, bm, mm.fb(), mapWeight, plane5AtomsESD));
	}

	if (not minimizer)
		throw runtime_error("Asymmetric unit with id " + asymID + " was not found");

//	float nudgeOffset = 0.05f;
//	if (vm.count("nudge-offset"))
//		nudgeOffset = vm["nudge-offset"].as<float>();
//	
//	minimizer->SetNudgeOffset(nudgeOffset);

	cout << "Initial score: " << minimizer->score() << endl
		 << "Initial rmsz values: " << endl;
	minimizer->printStats();

//	uint32_t iterations = 0;
//	if (vm.count("iterations"))
//		iterations = vm["iterations"].as<uint32_t>();

	minimizer->refine(true);
	
	cout << "Final score: " << minimizer->score() << endl
		 << "Final rmsz values: " << endl;
	minimizer->printStats();
	
	if (vm.count("stats"))
	{
		mmcif::StatsCollector collector(mm, structure, false /*electronScattering*/);
		auto r = collector.collect(asymID, resFirst, resLast);
	
		cout << "RESIDUE" << '\t'
			 << "RSR" << '\t'
			 << "SRSR" << '\t'
			 << "RSCCS" << '\t'
			 << "NGRID" << '\t'
			 << endl;
	
		for (auto i: r)
		{
			cout << fixed << setprecision(3)
				 << i.compID << '_' << i.asymID << '_' << i.seqID << '\t'
				 << i.RSR << '\t'
				 << i.SRSR << '\t'
				 << i.RSCCS << '\t'
				 << i.ngrid << '\t'
				 << endl;
		}
	}
		
	if (vm.count("xyzout"))
		f.save(vm["xyzout"].as<string>());
	
	return 0;
}
