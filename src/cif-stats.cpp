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
#include <regex>
#include <filesystem>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <gxrio.hpp>
#include <cif++.hpp>

namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;

class statsParser : public cif::sac_parser
{
  public:
	statsParser(std::istream& is)
		: sac_parser(is)
	{
	}
	
	void produce_datablock(const std::string& name) override
	{
	}
	
	void produce_category(const std::string& name) override
	{
	}
	
	void produce_row() override
	{
	}
	
	void produce_item(const std::string& category, const std::string& item, const std::string& value) override
	{
		size_t l = value.length();
		m_size_histogram[l] += 1;
	}
	
	std::map<size_t,size_t> m_size_histogram;
};

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-stats");
	visible_options.add_options()
		("input",	po::value<std::string>(),	"Input file")
		("help",								"Display help message")
		("version",								"Print version")
		("verbose,V",							"Verbose output");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",	po::value<int>(),			"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	
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
		exit(vm.count("help") ? 0 : 1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	gxrio::ifstream in(vm["input"].as<std::string>());
	if (not in.is_open())
		throw std::runtime_error("Could not open file " + vm["input"].as<std::string>());

	statsParser parser(in);
	parser.parse_file();

	size_t N = 0;
	for (const auto &[k, v] : parser.m_size_histogram)
	{
		N += v;
		std::cout << std::fixed << std::setw(6) << k << " : " << v << std::endl;
	}

	std::cout << std::endl
			  << "Total number of items: " << N << std::endl;

	std::cout << "sizeof a std::string is " << sizeof(std::string) << std::endl;

	return 0;
}

