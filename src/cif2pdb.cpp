/* include/cif++/Config.hpp.  Generated from Config.hpp.in by configure.  */
/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "cif-tools.hpp"

#include <sys/wait.h>

#include <fstream>
#include <chrono>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Cif2PDB.hpp"
#include "cif++/Structure.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

int pr_main(int argc, char* argv[])
{
	po::options_description desc("cif2pdb " + VERSION_STRING + " options");
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
		cout << argv[0] << " version " << VERSION_STRING << endl;
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
#warning "compile time PDB_DIR?"
	// if (not fs::exists(file) and regex_match(input, pdbIdRx))
	// 	file = fs::path(PDB_DIR) / "mmCIF" / input.substr(1, 2) / (input + ".cif.gz");
	
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

