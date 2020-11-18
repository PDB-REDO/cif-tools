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

#include <fstream>

#include <boost/program_options.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/CifUtils.hpp"

namespace po = boost::program_options;

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-validate [option...] file");
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("dict",	po::value<std::string>()->default_value("mmcif_pdbx_v50"),
												"The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file")
		("validate-links",						"Validate all links")
		("verbose,v",							"Verbose output");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input",	po::value<std::string>(),	"Input file")
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
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cerr << visible_options << std::endl;
		exit(vm.count("help") ? 0 : 1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	cif::File f;
	
	f.loadDictionary(vm["dict"].as<std::string>().c_str());
	
	if (vm.count("input") == 0)
		f.load(std::cin);
	else
		f.load(vm["input"].as<std::string>());

	int result = f.isValid() ? 0 : 1;

	if (vm.count("validate-links"))
		f.validateLinks();

	return result;
}

