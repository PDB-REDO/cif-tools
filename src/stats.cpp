/* 
   Created by: Maarten L. Hekkelman
   Date: woensdag 27 december, 2017
*/

#include "pdb-redo.h"

#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>


#include "zeep/el/element.hpp"

#include "cif++/BondMap.h"
#include "cif++/Statistics.h"

using namespace std;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace c = mmcif;

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("hklin",				po::value<string>(),	"mtz file")
		("recalc",										"Recalculate Fc from FP/SIGFP in mtz file")
		("aniso-scaling",		po::value<string>(),	"Anisotropic scaling (none/observed/calculated)")
		("no-bulk",										"No bulk correction")
		("xyzin",				po::value<string>(),	"coordinates file")
		("fomap",				po::value<string>(),	"Fo map file -- 2mFo - DFc")
		("dfmap",				po::value<string>(),	"difference map file -- 2(mFo - DFc)")
		("reshi",				po::value<float>(),		"High resolution")
		("reslo",				po::value<float>(),		"Low resolution")
		("sampling-rate",		po::value<float>(),		"Sampling rate")
		("electron-scattering",							"Use electron scattering factors")
		("no-edia",										"Skip EDIA score calculation")
		("output,o",			po::value<string>(),	"Write output to this file instead of stdout")
		("output-format",		po::value<string>(),	"Output format, can be either 'edstats' or 'json'")
		("use-auth-ids",								"Write auth_ identities instead of label_")
		("dict",				po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
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
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

//	po::variables_map vm;
//	po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
//	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		cerr << visible_options << endl;
		return 0;
	}
	
	if (vm.count("xyzin") == 0 or
		(vm.count("hklin") == 0 and (vm.count("fomap") == 0 or vm.count("dfmap") == 0)))
	{
		cerr << visible_options << endl;
		exit(1);
	}
	
	const set<string> kAnisoOptions{ "none", "calculated", "observed" };
	if (vm.count("aniso-scaling") and kAnisoOptions.count(vm["aniso-scaling"].as<string>()) == 0)
	{
		cerr << "Invalid option for aniso-scaling, allowed values are none, observed and calculated" << endl;
		exit(1);
	}
	
	if (vm.count("fomap") and (vm.count("reshi") == 0 or vm.count("reslo") == 0))
	{
		cerr << "The reshi and reslo parameters are required when using map files" << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// Load dict, if any
	
	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());

	mmcif::File f(vm["xyzin"].as<string>());
	mmcif::Structure structure(f);

	bool electronScattering = vm.count("electron-scattering") > 0;
	if (not electronScattering)
	{
		auto& exptl = f.data()["exptl"];
		electronScattering = not exptl.empty() and exptl.front()["method"] == "ELECTRON CRYSTALLOGRAPHY";
	}
	
	c::MapMaker<float> mm;
	
	if (vm.count("hklin"))
	{
		float samplingRate = 0.75;
		if (vm.count("sampling-rate"))
			samplingRate = vm["sampling-rate"].as<float>();
	
		if (vm.count("recalc"))
		{
			auto aniso = c::MapMaker<float>::as_None;
			if (vm.count("aniso-scaling"))
			{
				if (vm["aniso-scaling"].as<string>() == "observed")
					aniso = c::MapMaker<float>::as_Observed;
				else if (vm["aniso-scaling"].as<string>() == "calculated")
					aniso = c::MapMaker<float>::as_Calculated;
			}
			
			mm.calculate(
				vm["hklin"].as<string>(), structure, vm.count("no-bulk"), aniso, samplingRate, electronScattering);
		}
		else
			mm.loadMTZ(vm["hklin"].as<string>(), samplingRate);
	}
	else
	{
		float reshi = vm["reshi"].as<float>();
		float reslo = vm["reslo"].as<float>();
		
		mm.loadMaps(vm["fomap"].as<string>(), vm["dfmap"].as<string>(), reshi, reslo);
	}
	
	std::vector<mmcif::ResidueStatistics> r;
	
	if (vm.count("no-edia"))
	{
		mmcif::StatsCollector collector(mm, structure, electronScattering);
		r = collector.collect();
	}
	else
	{
		mmcif::BondMap bm(structure);

		mmcif::EDIAStatsCollector collector(mm, structure, electronScattering, bm);
		r = collector.collect();
	}

	ofstream of;
	if (vm.count("output"))
	{
		of.open(vm["output"].as<string>());
		if (not of.is_open())
		{
			cerr << "Could not open output file" << endl;
			exit(1);
		}
		cout.rdbuf(of.rdbuf());
	}

	if (vm.count("output-format") and vm["output-format"].as<string>() == "json")
	{
		using object = zeep::el::element;
		
		// vector<object> rs;
		
		// for (auto i: r)
		// {
		// 	object res;
		// 	res["asymID"] = i.asymID;
		// 	res["seqID"] = i.seqID;
		// 	res["compID"] = i.compID;
			
		// 	tuple<string,int,string,string> pdbID = structure.MapLabelToPDB(i.asymID, i.seqID, i.compID, i.authSeqID);
			
		// 	object pdb;
		// 	pdb["strandID"] = get<0>(pdbID);
		// 	pdb["seqNum"] = get<1>(pdbID);
		// 	pdb["compID"] = get<2>(pdbID);
		// 	pdb["insCode"] = get<3>(pdbID);
		// 	res["pdb"] = pdb;
			
		// 	res["RSR"] = i.RSR;
		// 	res["SRSR"] = i.SRSR;
		// 	res["RSCCS"] = i.RSCCS;
		// 	res["NGRID"] = i.ngrid;
		// 	res["EDIAm"] = i.EDIAm;
		// 	res["OPIA"] = i.OPIA;
			
		// 	rs.push_back(move(res));
		// }
		
		// object stats(rs);

		object stats;
		
		for (auto i: r)
		{
			tuple<string,int,string,string> pdbID = structure.MapLabelToPDB(i.asymID, i.seqID, i.compID, i.authSeqID);

			stats.emplace_back(object{
				{ "asymID", i.asymID },
				{ "seqID", i.seqID },
				{ "compID", i.compID },
				{
					"pdb", {
						{ "strandID", get<0>(pdbID) },
						{ "seqNum", get<1>(pdbID) },
						{ "compID", get<2>(pdbID) },
						{ "insCode", get<3>(pdbID) }
					}
				},
				{ "RSR", i.RSR },
				{ "SRSR", i.SRSR },
				{ "RSCCS", i.RSCCS },
				{ "NGRID", i.ngrid },
				{ "EDIAm", i.EDIAm },
				{ "OPIA", i.OPIA }
			});
		}
		
		cout << stats << endl;
	}
	else
	{
		cout << "RESIDUE" << '\t'
			 << "RSR" << '\t'
			 << "SRSR" << '\t'
			 << "RSCCS" << '\t'
			 << "NGRID" << '\t'
			 << "EDIAm" << '\t'
			 << "OPIA" << endl;
	
		bool writeAuth = vm.count("use-auth-ids");
	
		for (auto i: r)
		{
			string id;
			
			if (writeAuth)
			{
				tuple<string,int,string,string> pdbID = structure.MapLabelToPDB(i.asymID, i.seqID, i.compID, i.authSeqID);
				id = get<2>(pdbID) + '_' + get<0>(pdbID) + '_' + to_string(get<1>(pdbID)) + get<3>(pdbID);
			}
			else if (i.compID == "HOH")
				id = i.compID + '_' + i.asymID + '_' + i.authSeqID;
			else
				id = i.compID + '_' + i.asymID + '_' + to_string(i.seqID);
			
			cout << fixed << setprecision(3)
				 << id << '\t'
				 << i.RSR << '\t'
				 << i.SRSR << '\t'
				 << i.RSCCS << '\t'
				 << i.ngrid << '\t'
				 << i.EDIAm << '\t'
				 << setprecision(1) << i.OPIA << endl;
		}
	}
	
	return 0;
}
