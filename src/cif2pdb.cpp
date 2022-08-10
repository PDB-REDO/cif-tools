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
#include <gzstream/gzstream.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Cif2PDB.hpp"
#include "cif++/Structure.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif2pdb options input [output]");
	visible_options.add_options()
		("help,h",									"Display help message")
		("version",									"Print version")
		("verbose,v",								"Verbose output")
		("no-validate",								"Omit validation of the mmCIF file, forcing output in case of errors")
		("dict",		po::value<std::string>(),	"The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input",		po::value<std::string>(),	"Input file")
		("output,o",	po::value<std::string>(),	"Output file, default stdout")
		("debug,d",		po::value<int>(),			"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		write_version_string(std::cout, vm.count("verbose"));
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	std::string input = vm["input"].as<std::string>();
	std::regex pdbIdRx(R"(\d\w{3})");
	
	fs::path file = input;
	// #warning "compile time PDB_DIR?"
	// if (not fs::exists(file) and std::regex_match(input, pdbIdRx))
	// 	file = fs::path(PDB_DIR) / "mmCIF" / input.substr(1, 2) / (input + ".cif.gz");
	
	cif::File f;
	
	if (vm.count("dict"))
	{
		std::string dict = vm["dict"].as<std::string>();
		f.loadDictionary(dict.c_str());
	}
	else
		f.loadDictionary("mmcif_pdbx_v50");
	
	f.load(file);
	
	if (not vm.count("no-validate") and not f.isValid())
	{
		std::cerr << "This input mmCIF file is not valid";
		if (not cif::VERBOSE)
			std::cerr << ", use the --verbose option to see what errors were found" << std::endl;
		exit(1);
	}
	
	if (vm.count("output"))
	{
		file = vm["output"].as<std::string>();
		
		if (file.extension() == ".gz")
		{
			gzstream::ofstream out(file, std::ios_base::binary);
			WritePDBFile(out, f.front());
		}
		else
		{
			std::ofstream out(file.c_str(), std::ios_base::binary);
			WritePDBFile(out, f.front());
		}
	}
	else
		WritePDBFile(std::cout, f.front());
	
	return 0;	
}

