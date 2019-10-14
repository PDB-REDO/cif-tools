#include "libpr.h"

#include <fstream>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "cif++.h"
#include "cif-validator.h"
#include "pdb2cif.h"
#include "progress.h"
#include "cif-test.h"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

int drop(istream& is, set<string>& columns)
{
	cif::file in(is);
	
	for (auto c: columns)
	{
		string cat, item;
		tie(cat, item) = cif::split_tag_name(c);
		
		// loop over all datablocks
		for (auto& db: in)
		{
			auto& c = db[cat];
			if (not c.empty())
				c.drop(item);
		}
	}
	
	in.save(cout);
	
	return 0;
}

int cif_drop(int argc, char* argv[])
{
	po::options_description visible_options("cif-diff " + VERSION + " options file1 file2");
	visible_options.add_options()
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		("output,o",									"Write output to this file, default is to the terminal (stdout)")
		("column,c",	po::value<vector<string>>(),	"Column to drop, should be of the form '_category.item' with the leading underscore. Can be specified multiple times.");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",	po::value<string>(),		"Input file")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm.count("column") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	set<string> columns;
	for (auto cs: vm["column"].as<vector<string>>())
	{
		for (auto si = ba::make_split_iterator(cs, ba::token_finder(ba::is_any_of(",; "), ba::token_compress_on)); not si.eof(); ++si)
		{
			string c(si->begin(), si->end());
			ba::to_lower(c);
			columns.insert(c);
		}
	}

	if (cif::VERBOSE)
	{
		cerr << "Dropping columns:" << endl;
		for (auto c: columns)
			cerr << "    " << c << endl;
		cerr << endl;
	}
	
	fs::path file = vm["input"].as<string>();
	fs::ifstream is(file);
	if (not is.is_open())
	{
		cerr << "Could not open input file" << endl;
		exit(1);
	}

	ofstream f;
	if (vm.count("output"))
	{
		f.open(vm["output"].as<string>());
		if (not f.is_open())
		{
			cerr << "Could not open output file" << endl;
			exit(1);
		}
		cout.rdbuf(f.rdbuf());
	}
	
	return drop(is, columns);
}

