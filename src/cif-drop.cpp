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

#include <fstream>
#include <functional>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <cif++.hpp>

namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;

int drop(std::istream& is, std::set<std::string>& columns)
{
	cif::file in(is);
	
	for (auto c: columns)
	{
		std::string cat, item;
		std::tie(cat, item) = cif::split_tag_name(c);
		
		// loop over all datablocks
		for (auto& db: in)
		{
			auto& c = db[cat];
			if (not c.empty())
				c.drop(item);
		}
	}
	
	in.save(std::cout);
	
	return 0;
}

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-diff options file1 file2");
	visible_options.add_options()
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		("output,o",									"Write output to this file, default is to the terminal (stdout)")
		("column,c",	po::value<std::vector<std::string>>(),	"Column to drop, should be of the form '_category.item' with the leading underscore. Can be specified multiple times.");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",	po::value<std::string>(),		"Input file")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		write_version_string(std::cout, vm.count("verbose"));
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm.count("column") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	std::set<std::string> columns;
	for (auto cs: vm["column"].as<std::vector<std::string>>())
	{
		for (auto si = ba::make_split_iterator(cs, ba::token_finder(ba::is_any_of(",; "), ba::token_compress_on)); not si.eof(); ++si)
		{
			std::string c(si->begin(), si->end());
			ba::to_lower(c);
			columns.insert(c);
		}
	}

	if (cif::VERBOSE)
	{
		std::cerr << "Dropping columns:" << std::endl;
		for (auto c: columns)
			std::cerr << "    " << c << std::endl;
		std::cerr << std::endl;
	}
	
	fs::path file = vm["input"].as<std::string>();
	std::ifstream is(file);
	if (not is.is_open())
	{
		std::cerr << "Could not open input file" << std::endl;
		exit(1);
	}

	std::ofstream f;
	if (vm.count("output"))
	{
		f.open(vm["output"].as<std::string>());
		if (not f.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}
		std::cout.rdbuf(f.rdbuf());
	}
	
	return drop(is, columns);
}

