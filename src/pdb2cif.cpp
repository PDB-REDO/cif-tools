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

#include <gzstream/gzstream.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/PDB2Cif.hpp"
#include "cif++/Structure.hpp"
#include "cif++/Compound.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

int pr_main(int argc, char* argv[])
{
	std::string input;
	
	try
	{
		po::options_description visible_options("pdb2cif options input [output]");
		visible_options.add_options()
			("help,h",									"Display help message")
			("version",									"Print version")
			("verbose,v",								"Verbose output")
			("validate",								"Validate output file before writing")
			("dict",		po::value<std::string>(),	"Dictionary file containing restraints for residues in this specific target")
			;

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
		
		// Load dict, if any
		
		if (vm.count("dict"))
			c::CompoundFactory::instance().pushDictionary(vm["dict"].as<std::string>());
	
		input = vm["input"].as<std::string>();
		std::regex pdbIdRx(R"(\d\w{3})");
		
		fs::path file = input;

		std::unique_ptr<std::istream> in;

		if (file.extension() == ".gz")
		{
			gzstream::ifstream infile(file);

			if (not infile.is_open())
				throw std::runtime_error("Could not open file " + file.string());

			in.reset(new gzstream::ifstream(std::move(infile)));
		}
		else
		{
			std::ifstream infile(file);

			if (not infile.is_open())
				throw std::runtime_error("Could not open file " + file.string());

			in.reset(new std::ifstream(std::move(infile)));
		}
		
		cif::File f;
		ReadPDBFile(*in, f);
		
		if (vm.count("validate") and not f.isValid())
			throw std::runtime_error("The resulting mmCIF is not valid");
		
		if (vm.count("output"))
		{
			file = vm["output"].as<std::string>();

			if (file.extension() == ".gz")
			{
				gzstream::ofstream outfile(file, std::ios_base::binary);
				f.save(outfile);
			}
			else
			{
				std::ofstream outfile(file, std::ios_base::binary);
				f.save(outfile);
			}
		}
		else
			f.save(std::cout);
	}
	catch (const std::exception& ex)
	{
		if (not input.empty())
			std::cerr << "Error converting '" << input << '\'' << std::endl;
		throw;
	}

	return 0;	
}
