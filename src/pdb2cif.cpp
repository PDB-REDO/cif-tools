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
#include <stdexcept>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/PDB2Cif.hpp"
#include "cif++/Structure.hpp"
#include "cif++/Compound.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

int pr_main(int argc, char* argv[])
{
	string input;
	
	try
	{
		po::options_description desc("pdb2cif "s + VERSION_STRING + " options");
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
		
		// Load dict, if any
		
		if (vm.count("dict"))
			c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());
	
		input = vm["input"].as<string>();
		regex pdbIdRx(R"(\d\w{3})");
		
		fs::path file = input;
// #warning "compile time PDB_DIR?"
		// if (not fs::exists(file) and regex_match(input, pdbIdRx))
		// 	file = fs::path(PDB_DIR) / "pdb" / input.substr(1, 2) / ("pdb" + input + ".ent.gz");
		
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
