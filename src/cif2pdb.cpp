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
#include <filesystem>

#include <cfg.hpp>
#include <gxrio.hpp>

#include <cif++/pdb/io.hpp>

#include "revision.hpp"

namespace fs = std::filesystem;

int pr_main(int argc, char* argv[])
{
	auto &config = cfg::config::instance();

	config.init("usage: cif2pdb [options] inputfile [outputfile]",
		cfg::make_option("help,h",				"Display help message"),
		cfg::make_option("version",				"Print version"),
		cfg::make_option("verbose,v",			"Verbose output"),
		cfg::make_option("no-validate",			"Omit validation of the mmCIF file, forcing output in case of errors"),
		cfg::make_option<std::string>("dict",	"The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file"),
		cfg::make_hidden_option<int>("debug,d",	"Debug level (for even more verbose output)")
	);

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.operands().empty() or config.operands().size() > 2)
	{
		std::cerr << config << std::endl;
		exit(1);
	}

	cif::VERBOSE = config.has("verbose") != 0;
	if (config.has("debug"))
		cif::VERBOSE = config.get<int>("debug");
	
	// Load dict, if any
	
	if (config.has("dict"))
		cif::compound_factory::instance().push_dictionary(config.get<std::string>("dict"));

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	cif::VERBOSE = config.has("verbose") != 0;
	if (config.has("debug"))
		cif::VERBOSE = config.get<int>("debug");
	
	std::string input = config.operands().front();
	std::regex pdbIdRx(R"(\d\w{3})");
	
	fs::path file = input;
	// #warning "compile time PDB_DIR?"
	// if (not fs::exists(file) and std::regex_match(input, pdbIdRx))
	// 	file = fs::path(PDB_DIR) / "mmCIF" / input.substr(1, 2) / (input + ".cif.gz");
	
	gxrio::ifstream in(file);
	if (not in.is_open())
		throw std::runtime_error("Could not open file " + file.string());

	cif::file f(in);

	if (f.empty() or (not config.has("no-validate") and not f.is_valid()))
	{
		std::cerr << "This input mmCIF file is not valid";
		if (not cif::VERBOSE)
			std::cerr << ", use the --verbose option to see what errors were found" << std::endl;
		exit(1);
	}
	
	if (config.operands().size() == 2)
	{
		file = config.operands().back();

		gxrio::ofstream out(file);

		cif::pdb::write(out, f.front());
	}
	else
		cif::pdb::write(std::cout, f.front());
	
	return 0;	
}

