#include "pdb-redo.h"

#include <fstream>

#include <boost/program_options.hpp>

#include "cif++/mrsrc.h"

#include "cif++/Cif++.h"

using namespace std;
namespace po = boost::program_options;

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-validate [option...] file");
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("dict",	po::value<string>(),		"The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file")
		("validate-links",						"Validate all links")
		("verbose,v",							"Verbose output");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input",	po::value<string>(),		"Input file")
		("debug,d",	po::value<int>(),			"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		cerr << visible_options << endl;
		exit(vm.count("help") ? 0 : 1);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();

	cif::File f;
	
	if (vm.count("dict"))
	{
		string dict = vm["dict"].as<string>();
		f.loadDictionary(dict.c_str());
	}
	
	if (vm.count("input") == 0)
		f.load(cin);
	else
		f.load(vm["input"].as<string>());

	int result = f.isValid() ? 0 : 1;

	if (vm.count("validate-links"))
		f.validateLinks();

	return result;
}

