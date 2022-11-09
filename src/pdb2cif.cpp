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

#include <sys/wait.h>

#include <fstream>
#include <chrono>
#include <stdexcept>
#include <filesystem>

#include <cfg.hpp>

#include <cif++/pdb/io.hpp>
#include <cif++/gzio.hpp>

#include "revision.hpp"

namespace fs = std::filesystem;

int pr_main(int argc, char* argv[])
{
	std::string input;
	
	try
	{
		auto &config = cfg::config::instance();

		config.init("usage: pdb2cif [options] inputfile [outputfile]",
			cfg::make_option("help,h",				"Display help message"),
			cfg::make_option("version",				"Print version"),
			cfg::make_option("verbose,v",			"Verbose output"),
			cfg::make_option("validate",			"Validate output file before writing"),
			cfg::make_option<std::string>("dict",	"Dictionary file containing restraints for residues in this specific target"),
			cfg::make_hidden_option<int>("debug,d",	"Debug level (for even more verbose output)")
		);

		config.parse(argc, argv);
	
		if (config.has("version"))
		{
			write_version_string(std::cout, config.has("verbose"));
			exit(0);
		}

		if (config.has("help") or config.operands().empty() or config.operands().size() > 2)
		{
			std::cerr << config << std::endl;
			exit(config.has("help") ? 0 : 1);
		}
	
		cif::VERBOSE = config.has("verbose") != 0;
		if (config.has("debug"))
			cif::VERBOSE = config.get<int>("debug");
		
		// Load dict, if any
		
		if (config.has("dict"))
			cif::compound_factory::instance().push_dictionary(config.get<std::string>("dict"));
	
		input = config.operands().front();
		std::regex pdbIdRx(R"(\d\w{3})");
		
		fs::path file = input;

		cif::gzio::ifstream in(file);

		if (not in.is_open())
			throw std::runtime_error("Could not open file " + file.string());
		
		cif::file f = cif::pdb::read(in);
		
		if (config.has("validate") and not f.is_valid())
			throw std::runtime_error("The resulting mmCIF is not valid");
		
		if (config.operands().size() == 2)
		{
			file = config.operands().back();
			cif::gzio::ofstream out(file);
			f.save(out);
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
