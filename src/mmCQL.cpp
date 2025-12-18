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

#include "cif++/datablock.hpp"
#include "revision.hpp"

#include <cif++.hpp>
#include <cif++/cql.hpp>
#include <cif++/gzio.hpp>
#include <exception>
#include <fstream>
#include <linenoise.h>
#include <mcfp/mcfp.hpp>

// -----------------------------------------------------------------------

int pr_main(int argc, char *argv[])
{
	auto &config = mcfp::config::instance();

	config.init("mmCQL [options] input [output]",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,V", "Verbose output, repeat to increase verbosity level"),

		mcfp::make_option("force,F", "Force writing of output file, even if it is the same as the input file"),

		mcfp::make_option<std::string>("file,f", "Read SQL commands from file"),
		mcfp::make_option<std::string>("command,c", "Single SQL script to execute"),

		mcfp::make_option<std::string>("data-block,D", "Datablock to use, default is first"));

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("file") and config.has("command"))
	{
		std::cerr << "Please specify either 'file' or 'command', not both\n";
		exit(1);
	}

	if (config.has("help") or config.operands().empty() or config.operands().size() > 2)
	{
		std::cerr << config << '\n';
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");

	if (config.operands().size() == 2 and config.operands().front() == config.operands().back() and not config.has("force"))
	{
		std::cerr << "Cowardly refusing to overwrite input file (specify --force to force overwriting)\n";
		exit(1);
	}

	cif::gzio::ifstream in(config.operands().front());
	if (not in.is_open())
		throw std::runtime_error("Could not open file " + config.operands().front());

	cif::file file{ in };

	cif::datablock &db = config.has("data-block") ? file[config.get("data-block")] : file.front();
	db.load_dictionary();

	cif::cql::connection connection(db);

	if (config.has("file"))
	{
		std::ifstream cmdFile(config.get<std::string>("file"));
		if (not cmdFile.is_open())
			throw std::runtime_error("Failed to open command file " + config.get<std::string>("file"));

		cif::cql::transaction tx(connection);

		std::stringstream ss;
		ss << cmdFile.rdbuf();
		auto r = tx.exec(ss.str());
		if (not r.empty())
			std::cout << r << "\n";
		tx.commit();
	}
	else if (config.has("command"))
	{
		cif::cql::transaction tx(connection);
		auto r = tx.exec(config.get("command"));
		if (not r.empty())
			std::cout << r << "\n";
		tx.commit();
	}
	else
	{
		const char *line;

		linenoiseHistoryLoad(".mmcql-history");

		while ((line = linenoise("cql> ")) != nullptr)
		{
			try
			{
				cif::cql::transaction tx(connection);
				auto r = tx.exec(line);
				std::cout << r << "\n";
				tx.commit();
			}
			catch (const std::exception &ex)
			{
				std::cout << "Error executing statement(s): " << ex.what() << "\n";
			}

			linenoiseHistoryAdd(line);

			free((void *)line);
		}

		linenoiseHistorySave(".mmcql-history");
	}

	if (config.operands().size() == 2)
	{
		cif::gzio::ofstream out(config.operands().back());
		if (not out.is_open())
			throw std::runtime_error("Could not open output file " + config.operands().back());

		file.save(out);
	}

	return 0;
}
