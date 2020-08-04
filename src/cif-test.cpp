#include "libpr.h"

#include <sys/wait.h>

#include <fstream>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include "cif++.h"
#include "pdb2cif.h"
#include "cif2pdb.h"
#include "cif-test.h"
#include "peptidedb.h"

#include "mmcif/structure.h"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

int cif_test(int argc, char* argv[])
{
	po::options_description visible_options("cif-test " + VERSION_STRING + " options file1 [file2...]" );
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("dict",		po::value<string>(),	"Load dictionary, default is mmcif_ddl. If this value is a local file, the content of the file is used, otherwise it uses a compiled in resource. Currently mmcif_ddl (version 2.1.6) and mmcif_pdbx (version 5.286) are supported.")
		("validate",							"Validate the content of the file using the default or specified dictionary")
		("verbose,v",							"Verbose output")
		("print,p",								"Dump internal CIF data");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",		po::value<vector<string>>(),	"Input files")
		("width,w",		po::value<int>(),				"Width for text wrap")
		("debug,d",		po::value<int>(),				"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", -1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("width"))
	{
		int w = vm["width"].as<int>();

const char s[] = R"(1. HYDROGENS HAVE BEEN ADDED IN THE RIDING POSITIONS.
2. A MET-INHIBITION PROTOCOL WAS USED FOR SELENOMETHIONINE INCORPORATION DURING PROTEIN EXPRESSION. THE OCCUPANCY OF THE SE ATOMS IN THE MSE RESIDUES WAS REDUCED TO 0.75 FOR THE REDUCED SCATTERING POWER DUE TO PARTIAL S-MET INCORPORATION.
3. ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY.
4. THERE IS UNMODELED DENSITY NEAR SER 10 IN EACH CHAIN, WHICH IS THE PUTATIVE GLUTATHIONE BINDING SITE.
5. PO4 AND EDO WERE MODELED BASED ON CRYSTALLIZATION AND CRYOPROTECTION CONDITIONS.
)";
		
		auto lines = cif::word_wrap(s, w);
		
		for (auto line: lines)
			cout << "  | " << line << string(w - line.length(), ' ') << " |" << endl;
				
		exit(0);
	}

	if (vm.count("version"))
	{
		cout << argv[0] << " version " VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	for (fs::path input: vm["input"].as<vector<string>>())
	{
		cout << "validating " << input << endl;
		
		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("No such file");
		
		io::filtering_stream<io::input> in;
		
		if (input.extension() == ".bz2")
		{
			in.push(io::bzip2_decompressor());
			input = input.stem();
		}
		else if (input.extension() == ".gz")
		{
			in.push(io::gzip_decompressor());
			input = input.stem();
		}
		
		in.push(infile);
		
//		using namespace std::chrono;
//		
//		auto start = steady_clock::now();
//		cout << "start reading" << endl;
	
		cif::file file(in);
	
//		auto end = steady_clock::now();
//		cout << "reading done" << endl;
//		
//		duration<double> d = duration_cast<duration<double>>(end - start);
//		
//		cout << "Time elapsed in seconds: " << d.count() << endl;
		
		if (vm.count("validate"))
		{
			if (vm.count("dict"))
			{
				string dict = vm["dict"].as<string>();
				if (fs::exists(dict))
				{
					ifstream dictFile(dict);
					file.load_dictionary(dictFile);
				}
				else
					file.load_dictionary(dict.c_str());
			}
			else
				file.load_dictionary();
		
			file.validate();
		}

		if (vm.count("print"))
			file.save(cout);
	}
//	
//	auto& e = file.first_datablock();
////	cout << "ID of datablock: " << e.first_item("_datablock.id").as<string>() << endl;
//
//	cout << "ID of datablock: " << e["entry"].get_first_item("id") << endl;
//	
//	for (auto& cat: e)
//		cout << "category: " << cat.name() << endl;
//	
//	auto& atom_site = e["atom_site"];
//	for (auto& atom: atom_site)
//	{
////		if (atom["group_PDB"].as<string>() != "HETATM")
////			continue;
//		
//		if (atom["label_comp_id"] != "HOH")
//			continue;
//		
//		float x, y, z;
//		cif::tie(x, y, z) = atom.get("Cartn_x", "Cartn_y", "Cartn_z");
//		
//		cout << "water at: " << x << ',' << y << ',' << z << endl;
//	}
	
	return 0;
}

int pdb2cif_test(int argc, char* argv[])
{
	po::options_description desc("pdb2cif-test " + VERSION_STRING + " options");
	desc.add_options()
		("input,i",		po::value<string>(),	"Input file")
		("output,o",	po::value<string>(),	"Output file, default stdout")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)")
		("dict",		po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
		;

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " VERSION_STRING << endl;
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
		CompoundFactory::instance().PushDictionary(vm["dict"].as<string>());

	fs::path file = vm["input"].as<string>();
	
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

	cif::file f;
	ReadPDBFile(in, f);
	
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
	
	return 0;	
}

int cif2pdb_test(int argc, char* argv[])
{
	po::options_description desc("cif2pdb-test " + VERSION_STRING + " options");
	desc.add_options()
		("input,i",		po::value<string>(),	"Input file")
		("output,o",	po::value<string>(),	"Output file, default stdout")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " VERSION_STRING << endl;
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
	
	fs::path file = vm["input"].as<string>();
	
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

	cif::file f(in);
	
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

int pdb_diff_test(int argc, char* argv[])
{
	po::options_description visible_options("cif2pdb-diff " + VERSION_STRING + " [options] PDB-ID");
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("input,i",		po::value<string>(),	"PDB ID")
		("data_dir",	po::value<string>(),	"Directory containing both mmCIF and pdb directories");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

#warning "compile time PDB_DIR?"
	fs::path baseDir;// = PDB_DIR;
	if (vm.count("data_dir"))
		baseDir = vm["data_dir"].as<string>();
	if (not fs::exists(baseDir))
		throw runtime_error("base dir " + baseDir.string() + " does not exist");
	
	string pdbID = vm["input"].as<string>();
	ba::to_lower(pdbID);
	if (regex_match(pdbID, regex(R"([0-9][0-9a-z]{3})")) == false)
		throw runtime_error("Invalid pdb ID");
	
	fs::path mmCIFFile = baseDir / "mmCIF" / pdbID.substr(1, 2) / (pdbID + ".cif.gz");
	if (not fs::exists(mmCIFFile))
		throw runtime_error("mmCIF file " + mmCIFFile.string() + " does not exist");

	fs::path pdbFile = baseDir / "pdb" / pdbID.substr(1, 2) / ("pdb" + pdbID + ".ent.gz");
	if (not fs::exists(pdbFile))
		throw runtime_error("PDB file " + pdbFile.string() + " does not exist");
	
	ifstream infile(pdbFile.c_str(), ios_base::in | ios_base::binary);
	if (not infile.is_open())
		throw runtime_error("Could not open file " + pdbFile.string());

	io::filtering_stream<io::input> in;
	if (pdbFile.extension() == ".bz2")
		in.push(io::bzip2_decompressor());
	else if (pdbFile.extension() == ".gz")
		in.push(io::gzip_decompressor());
	in.push(infile);

	// temp files for vimdiff
	// first, the original file
	
	char generated[] = "/tmp/cif2pdb-diff-generated-XXXXXX.pdb", original[] = "/tmp/cif2pdb-diff-original-XXXXXX.pdb";
	int fd[2];
	
	if ((fd[0] = mkstemps(generated, 4)) < 0 or (fd[1] = mkstemps(original, 4)) < 0)
	{
		cerr << "Error creating temp files:  " << strerror(errno) << endl;
		exit(1);
	}

	io::file_descriptor_sink gen(fd[0], io::close_handle), orig(fd[1], io::close_handle);

	io::filtering_stream<io::output> out;
	out.push(orig);
	io::copy(in, out);

	// Next the converted cif file

	ifstream infile2(mmCIFFile.c_str(), ios_base::in | ios_base::binary);
	if (not infile2.is_open())
		throw runtime_error("Could not open file " + mmCIFFile.string());
	
	io::filtering_stream<io::input> in2;
	if (mmCIFFile.extension() == ".bz2")
		in2.push(io::bzip2_decompressor());
	else if (mmCIFFile.extension() == ".gz")
		in2.push(io::gzip_decompressor());
	in2.push(infile2);

	cif::file cf(in2);

	io::filtering_stream<io::output> out2;
	out2.push(gen);
	
	WritePDBFile(out2, cf);
	
	const char* const n_argv[] = {
		"/usr/bin/vimdiff",
		"-c", "set diffopt-=iwhite",
		original,
		generated,
		nullptr
	};
	
	int pid = fork();
	
	if (pid <= 0)
	{
		if (execv(n_argv[0], const_cast<char*const*>(n_argv)) < 0)
			cerr << "Failed to execute vimdiff" << endl;
		exit(1);
	}
	
	int status;
	waitpid(pid, &status, 0);
	
	if (WIFEXITED(status))
	{
		unlink(generated);
		unlink(original);
	}	

	return 0;	
}

int pdb2pdb_diff_test(int argc, char* argv[])
{
	po::options_description visible_options("cif2pdb-diff " + VERSION_STRING + " [options] PDB-ID");
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("input,i",		po::value<string>(),	"PDB ID or file")
		("icase",								"Ignore case (vimdiff option)")
		("iwhite",								"Ignore whitespace (vimdiff option)")
		("data_dir",	po::value<string>(),	"Directory containing both mmCIF and pdb directories");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

#warning "compile time PDB_DIR?"
	fs::path baseDir;// = PDB_DIR;
	if (vm.count("data_dir"))
		baseDir = vm["data_dir"].as<string>();
	if (not fs::exists(baseDir))
		throw runtime_error("base dir " + baseDir.string() + " does not exist");
	
	fs::path pdbFile;
	
	string pdbID = vm["input"].as<string>();
	if (regex_match(pdbID, regex(R"([0-9][0-9a-z]{3})", regex_constants::icase)))
	{
		ba::to_lower(pdbID);
	
		pdbFile = baseDir / "pdb" / pdbID.substr(1, 2) / ("pdb" + pdbID + ".ent.gz");
	}
	else
		pdbFile = pdbID;

	if (not fs::exists(pdbFile))
		throw runtime_error("PDB file " + pdbFile.string() + " does not exist");
	
	ifstream infile(pdbFile.c_str(), ios_base::in | ios_base::binary);
	if (not infile.is_open())
		throw runtime_error("Could not open file " + pdbFile.string());

	io::filtering_stream<io::input> in;
	if (pdbFile.extension() == ".bz2")
		in.push(io::bzip2_decompressor());
	else if (pdbFile.extension() == ".gz")
		in.push(io::gzip_decompressor());
	in.push(infile);

	// temp files for vimdiff
	// first, the original file
	
	char generated[] = "/tmp/cif2pdb-diff-generated-XXXXXX.pdb", original[] = "/tmp/cif2pdb-diff-original-XXXXXX.pdb";
	int fd[2];
	
	if ((fd[0] = mkstemps(generated, 4)) < 0 or (fd[1] = mkstemps(original, 4)) < 0)
	{
		cerr << "Error creating temp files:  " << strerror(errno) << endl;
		exit(1);
	}

	io::file_descriptor_sink gen(fd[0], io::close_handle), orig(fd[1], io::close_handle);

	io::filtering_stream<io::output> out;
	out.push(orig);
	io::copy(in, out);

	// Next the converted cif file

	ifstream infile2(pdbFile.c_str(), ios_base::in | ios_base::binary);
	if (not infile2.is_open())
		throw runtime_error("Could not open file " + pdbFile.string());
	
	io::filtering_stream<io::input> in2;
	if (pdbFile.extension() == ".bz2")
		in2.push(io::bzip2_decompressor());
	else if (pdbFile.extension() == ".gz")
		in2.push(io::gzip_decompressor());
	in2.push(infile2);

	cif::file cf;
	ReadPDBFile(in2, cf);

	io::filtering_stream<io::output> out2;
	out2.push(gen);
	
	WritePDBFile(out2, cf);
	
	vector<const char*> n_argv = {
		"/usr/bin/vimdiff"
	};
	
	if (vm.count("icase"))
	{
		n_argv.push_back("-c");
		n_argv.push_back("set diffopt+=icase");
	}
	
	if (vm.count("iwhite"))
	{
		n_argv.push_back("-c");
		n_argv.push_back("set diffopt-=iwhite");
	}

	n_argv.push_back(original);
	n_argv.push_back(generated);
	n_argv.push_back(nullptr);
	
	int pid = fork();
	
	if (pid <= 0)
	{
		if (execv(n_argv[0], const_cast<char*const*>(n_argv.data())) < 0)
			cerr << "Failed to execute vimdiff" << endl;
		exit(1);
	}
	
	int status;
	waitpid(pid, &status, 0);
	
	if (WIFEXITED(status))
	{
		unlink(generated);
		unlink(original);
	}	

	return 0;	
}
