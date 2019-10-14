#include "pdb-redo.h"

#include <sys/wait.h>

#include <fstream>
#include <chrono>
#include <stdexcept>

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "cif++/Cif++.h"
#include "cif++/PDB2Cif.h"
#include "cif++/Structure.h"
#include "cif++/Compound.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

int pr_main(int argc, char* argv[])
{
	string input;
	
	try
	{
		po::options_description desc("pdb2cif " + VERSION + " options");
		desc.add_options()
			("input,i",		po::value<string>(),	"Input file")
			("output,o",	po::value<string>(),	"Output file, default stdout")
			("help,h",								"Display help message")
			("version",								"Print version")
			("verbose,v",							"Verbose output")
			("validate",							"Validate output file before writing")
			("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)")
			("dict",		po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
			;
	
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
		
		// Load dict, if any
		
		if (vm.count("dict"))
			c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());
	
		input = vm["input"].as<string>();
		regex pdbIdRx(R"(\d\w{3})");
		
		fs::path file = input;
		if (not fs::exists(file) and regex_match(input, pdbIdRx))
			file = fs::path(PDB_DIR) / "pdb" / input.substr(1, 2) / ("pdb" + input + ".ent.gz");
		
		ifstream infile(file.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("Could not open file " + file.string());
	
		io::filtering_stream<io::input> in;
	
		if (file.extension() == ".bz2")
		{
			in.push(io::bzip2_decompressor());
			file = file.stem();
		}
		else if (file.extension() == ".gz")
		{
			in.push(io::gzip_decompressor());
			file = file.stem();
		}
		
		in.push(infile);
	
		cif::File f;
		ReadPDBFile(in, f);
		
		if (vm.count("validate") and not f.isValid())
			throw runtime_error("The resulting mmCIF is not valid");
		
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
			
			f.save(out);
		}
		else
			f.save(cout);
	}
	catch (const exception& ex)
	{
		if (not input.empty())
			cerr << "Error converting '" << input << '\'' << endl;
		throw;
	}

	return 0;	
}
