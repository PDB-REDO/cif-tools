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

#include <fstream>

#include <cfg.hpp>
#include <cif++.hpp>

#include "revision.hpp"

int pr_main(int argc, char *argv[])
{
	auto &config = cfg::config::instance();

	config.init(
		cfg::make_option("help,h", "Display help message"),
		cfg::make_option("version", "Print version"),
		cfg::make_option<std::string>("dict", "mmcif_pdbx.dic", "The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file"),
		cfg::make_option("validate-links", "Validate all links"),
		cfg::make_option("verbose,v", "Verbose output"),
		cfg::make_option<int>("debug,d", "Debug level (for even more verbose output)"));

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().empty())
	{
		std::cerr << "cif-validate [options...] file" << std::endl
				  << std::endl
				  << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");
	if (config.has("debug"))
		cif::VERBOSE = config.get<int>("debug");

	cif::file f(config.operands().front());

	if (config.count("dict"))
		f.load_dictionary(config.get<std::string>("dict"));

	int result = f.is_valid() ? 0 : 1;

	if (config.has("validate-links"))
		f.validate_links();

	return result;
}
