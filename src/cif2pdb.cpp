#include "pdb-redo.h"

#include <sys/wait.h>

#include <fstream>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "cif++/Cif++.h"
#include "cif++/Cif2PDB.h"
#include "cif++/Structure.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

int pr_main(int argc, char* argv[])
{
	po::options_description desc("cif2pdb " + VERSION + " options");
	desc.add_options()
		("input,i",		po::value<string>(),	"Input file")
		("output,o",	po::value<string>(),	"Output file, default stdout")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("no-validate",							"Omit validation of the mmCIF file, forcing output in case of errors")
		("dict",		po::value<string>(),	"The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << desc << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	string input = vm["input"].as<string>();
	regex pdbIdRx(R"(\d\w{3})");
	
	fs::path file = input;
	if (not fs::exists(file) and regex_match(input, pdbIdRx))
		file = fs::path(PDB_DIR) / "mmCIF" / input.substr(1, 2) / (input + ".cif.gz");
	
	cif::File f;
	
	if (vm.count("dict"))
	{
		string dict = vm["dict"].as<string>();
		f.loadDictionary(dict.c_str());
	}
	else
		f.loadDictionary("mmcif_pdbx");
	
	f.load(file);
	
	if (not vm.count("no-validate") and not f.isValid())
	{
		cerr << "This input mmCIF file is not valid";
		if (not cif::VERBOSE)
			cerr << ", use the --verbose option to see what errors were found" << endl;
		exit(1);
	}
	
	if (vm.count("output"))
	{
		file = vm["output"].as<string>();
		
		ofstream outfile(file.c_str(), ios_base::out | ios_base::binary);
		io::filtering_stream<io::output> out;
		
		if (file.extension() == ".gz")
			out.push(io::gzip_compressor());
		else if (file.extension() == ".bz2")
			out.push(io::bzip2_compressor());
		
		out.push(outfile);
		
		WritePDBFile(out, f);
	}
	else
		WritePDBFile(cout, f);
	
	return 0;	
}

