/* 
   Created by: Maarten L. Hekkelman
   Date: woensdag 07 maart, 2018
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
#include "cif++/CifUtils.h"

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

vector<string> CollectCandidateIds(bool includeSugars)
{
	const char* clibd_mon = getenv("CLIBD_MON");
	if (clibd_mon == nullptr)
		throw runtime_error("CLIBD_MON is not defined, please source the CCP4 environment");

	fs::path monomersDir(clibd_mon);
	
	vector<string> ids;
	
	for (fs::directory_iterator mp(monomersDir); mp != fs::directory_iterator(); ++mp)
	{
		if (not is_directory(*mp))
			continue;
			
		for (fs::directory_iterator mps(mp->path()); mps != fs::directory_iterator(); ++mps)
		{
			if (is_directory(*mps))
				continue;
			
			string n = mps->path().filename().string();
			if (ba::ends_with(n, ".cif"))
			{
				n = n.substr(0, n.length() - 4);
				if (n == "COM_COM" or n == "CON_CON" or n == "PRN_PRN")
					n = n.substr(0, 3);
				ids.push_back(n);
			}
		}
	}
	
	sort(ids.begin(), ids.end());
	
	cif::Progress p1(ids.size(), "init");
//	for (auto i = ids.begin(); i != ids.end(); ++i)
//	{
//		p1.message(*i);
//		
//		try
//		{
//			auto l = mmcif::Compound::create(*i);
//			if (l == nullptr or l->chiralCentres().empty())
//				i = ids.erase(i) - 1;
//		}
//		catch (const exception& ex)
//		{
//			cerr << ex.what() << endl;
//			i = ids.erase(i) - 1;
//		}
//
//		p1.consumed(1);
//	}

	boost::thread_group t;
	size_t N = boost::thread::hardware_concurrency();
	atomic<size_t> next(0);
	
	boost::mutex m;
	boost::barrier bar(N);

	for (size_t i = 0; i < N; ++i)
	{
		t.create_thread([&]()
		{
			set<string> failed;
			
			for (;;)
			{
				size_t i = next++;
				
				if (i >= ids.size())
					break;
				
				try
				{
					
					auto l = mmcif::Compound::create(ids[i]);
					if (l == nullptr or
						(l->isSugar() and not includeSugars) or
						l->chiralCentres().empty())
					{
						failed.insert(ids[i]);
					}
				}
				catch (const exception& ex)
				{
					cerr << ex.what() << endl;
					failed.insert(ids[i]);
				}
		
				p1.consumed(1);
			}
			
			bar.wait();
			
			boost::mutex::scoped_lock l(m);
			
			vector<string> nIds;
			nIds.reserve(ids.size() - failed.size());

			set_difference(ids.begin(), ids.end(), failed.begin(), failed.end(),
				back_inserter(nIds));
			
			swap(ids, nIds);
		});
	}

	t.join_all();
	
	return ids;
}

// --------------------------------------------------------------------

void CompareCompounds(const vector<string>& ids, ostream& report)
{
	boost::thread_group t;
	size_t N = boost::thread::hardware_concurrency();
	atomic<size_t> next(0);
	
	boost::mutex m;
	
	vector<tuple<int,int>> idSets;

	for (size_t i = 0; i + 1 < ids.size(); ++i)
	{
		for (size_t j = i + 1; j < ids.size(); ++j)
			idSets.emplace_back(i, j);
	}

	cif::Progress p2(idSets.size(), "compare");

	for (size_t i = 0; i < N; ++i)
	{
		t.create_thread([&]()
		{
			set<string> failed;
			
			for (;;)
			{
				size_t i = next++;
				
				if (i >= idSets.size())
					break;
				
				int a, b;
				tie(a, b) = idSets[i];
				
				auto* ca = mmcif::Compound::create(ids[a]);
				auto* cb = mmcif::Compound::create(ids[b]);
				
				p2.consumed(1);

				if (ca->isIsomerOf(*cb))
				{
					boost::mutex::scoped_lock l(m);
					report << ids[a] << " is isomer of " << ids[b] << endl;
				}
			}
		});
	}

	t.join_all();
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("help,h",										"Display help message")
		("version",										"Print version")
		("with-sugar",									"Include sugars in isomer sets")
		("output,o", po::value<string>(),				"Output file")
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
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
	
	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();

#if 1
	auto ids = CollectCandidateIds(vm.count("with-sugar") > 0);
	
	cout << "Comparing " << ids.size() << " entries all against all" << endl;
	
	if (vm.count("output"))
	{
		ofstream f(vm["output"].as<string>(), ios_base::out);
		if (not f.is_open())
			throw runtime_error("Could not open output file " + vm["output"].as<string>());
		
		CompareCompounds(ids, f);
	}
	else
		CompareCompounds(ids, cout);
#else
	CompareCompounds({ "OKA", "XT2" }, cout);
#endif
	return 0;
}
