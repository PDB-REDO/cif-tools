/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "mrsrc.hpp"
#include "revision.hpp"

#include <cerrno>
#include <cif++.hpp>
#include <cif++/category.hpp>
#include <cif++/cql.hpp>
#include <cif++/datablock.hpp>
#include <cif++/gzio.hpp>
#include <cif++/text.hpp>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <istream>
#include <mcfp/mcfp.hpp>
#include <memory>
#include <ranges>
#include <readln.hpp>
#include <sstream>
#include <string_view>
#include <sys/poll.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <system_error>
#include <unistd.h>
#include <utility>

// --------------------------------------------------------------------
// Globals to store settings

bool gOutputModeOnlyRows = false;
cif::category::output_format gOutputFormat = cif::category::output_format::column;
std::ofstream gOutputFile;

// --------------------------------------------------------------------

uint32_t get_terminal_height()
{
	uint32_t result = 24;

	if (::isatty(STDOUT_FILENO))
	{
		struct winsize w;
		::ioctl(0, TIOCGWINSZ, &w);
		result = w.ws_row;
	}

	return result;
}

// --------------------------------------------------------------------
// Data formatting

void writeResult(cif::category &cat, std::ostream &os)
{
	cat.write(os, gOutputFormat, {}, false);
}

void writeResult(cif::category &cat)
{
	if (gOutputFile.is_open())
		writeResult(cat, gOutputFile);
	else
		writeResult(cat, std::cout);
}

// std::istream *writeResult(cif::category &cat, OutputFormat format)
// {
// 	auto result = std::make_unique<std::stringstream>();

// 	if (format == OutputFormat::CIF)
// 		cat.write(*result, {}, false);
// 	else
// 	{

// 	}

// 	return result.release();
// }

// --------------------------------------------------------------------

void showPagerForData(std::istream &is)
{
	std::error_code ec;
	auto tmpfile = std::filesystem::temp_directory_path(ec) / std::format("mmcql-fifo-{}", getpid());
	if (ec == std::errc{})
	{
		int r = mkfifo(tmpfile.c_str(), 0600);
		if (r != 0)
			ec = std::error_code(errno, std::system_category());
	}

	if (ec)
		std::cout << is.rdbuf();
	else
	{
		switch (auto pid = fork())
		{
			case -1:
				std::cerr << "fork failed: " << std::error_code(errno, std::system_category()).message() << "\n";
				break;
			case 0:
				execlp("/bin/sh", "sh", "-c", "pager -d -eF /tmp/mmcql-fifo", nullptr);
				std::cerr << "exec of pager failed: " << std::error_code(errno, std::system_category()).message() << "\n";
				exit(-1);
				break;
			default:
				std::ofstream of("/tmp/mmcql-fifo");
				of << is.rdbuf();
				of.close();
				waitpid(pid, nullptr, 0);
				std::filesystem::remove(tmpfile, ec);
				break;
		}
	}
}

// --------------------------------------------------------------------

void displayHelp()
{
	std::cout << R"(You are using mmcql, a command-line interface for CIF files using SQLite as engine.
Type:  \copyright for distribution terms
       \h for help with SQL commands
       \? for help with mmcql commands
       \q to quit
)";
	//    \g or terminate with semicolon to execute query
}

void displayOptions();

// --------------------------------------------------------------------

enum class CommandCategory
{
	General,
	Help,
	Input_Output,
	Informational,
	Formatting,
	Datablock
};

std::map<CommandCategory, std::string> kCommandCategoryLabels{
	{ CommandCategory::General, "General" },
	{ CommandCategory::Help, "Help" },
	{ CommandCategory::Input_Output, "Input/Output" },
	{ CommandCategory::Informational, "Informational" },
	{ CommandCategory::Formatting, "Formatting" },
	{ CommandCategory::Datablock, "Datablock" }
};

struct BackslashCommand
{
	CommandCategory mCategory;
	std::string_view mCommand;
	std::string_view mDesc;
	std::function<void(std::string_view arg)> mFunc;
};

std::vector<BackslashCommand> gBackslashCommands{
	{ CommandCategory::General,
		"\\copyright",
		"show copyright and usage", [](std::string_view)
		{
			mrsrc::istream license("LICENSE");
			showPagerForData(license);
		} },
	{ CommandCategory::General, "\\?", "show help on backslash commands", [](std::string_view)
		{
			displayOptions();
		} },
	{ CommandCategory::General, "\\h", "show help on SQL syntax", [](std::string_view)
		{
			std::cout << "This will eventually be the SQL help page\n";
		} },
	{ CommandCategory::Formatting, "\\format", "Select one of: cif, csv, tsv, list, column, markdown, table", [](std::string_view fmt)
		{
			std::string f{ fmt };
			cif::trim(f);

			if (f == "cif")
				gOutputFormat = cif::category::output_format::cif;
			else if (f == "csv")
				gOutputFormat = cif::category::output_format::csv;
			else if (f == "tsv")
				gOutputFormat = cif::category::output_format::tsv;
			else if (f == "list")
				gOutputFormat = cif::category::output_format::list;
			else if (f == "column")
				gOutputFormat = cif::category::output_format::column;
			else if (f == "markdown")
				gOutputFormat = cif::category::output_format::markdown;
			else if (f == "table")
				gOutputFormat = cif::category::output_format::table;
			else
			 	std::cout << "Unknown output format: " << f << "\n";
		} }
};

void displayOptions()
{
	std::stringstream ss;
	for (const auto &[cat, label] : kCommandCategoryLabels)
	{
		for (bool first = true; auto cmd : gBackslashCommands | std::views::filter([cat](BackslashCommand &cmd)
																	{ return cmd.mCategory == cat; }))
		{
			if (std::exchange(first, false))
				ss << label << "\n\n";
			ss << std::format("  {:10}  {}\n", cmd.mCommand, cmd.mDesc);
		}
	}
	showPagerForData(ss);
}

// void completion(const char *buf, linenoiseCompletions *lc)
// {
// 	using namespace std::literals;

// 	std::string_view bufsv(buf);

// 	for (auto bc : gBackslashCommands)
// 	{
// 		for (auto l = bc.mCommand.size(); l > 0; --l)
// 		{
// 			if (bufsv.ends_with(bc.mCommand.substr(0, l)))
// 			{
// 				linenoiseAddCompletion(lc, bc.mCommand.data());
// 				break;
// 			}
// 		}
// 	}
// }

// --------------------------------------------------------------------

class MMCQLApplication
{
  public:
	MMCQLApplication();

	void loop();
	void loadCifFile(std::string_view f);
	void loadDatablock(std::string_view d);

	~MMCQLApplication();

  private:
	void showPagerForData(cif::category &cat);

	void showCategories(std::string_view);
	void showDatablocks(std::string_view);

	void processCommand(std::string_view cmd, std::string_view args);

	std::filesystem::path m_filename;
	std::unique_ptr<cif::file> m_file;
	std::string m_dbname;
	std::unique_ptr<cif::cql::connection> m_connection;
	bool m_modified = false;
};

// --------------------------------------------------------------------

MMCQLApplication::MMCQLApplication()
{
	gBackslashCommands.emplace_back(
		CommandCategory::Input_Output,
		"\\file", "load new mmCIF file", [this](std::string_view f)
		{ this->loadCifFile(f); });
	gBackslashCommands.emplace_back(
		CommandCategory::Datablock,
		"\\dd", "show datablocks in file", [this](std::string_view a)
		{ this->showDatablocks(a); });
	gBackslashCommands.emplace_back(
		CommandCategory::Datablock,
		"\\dc", "show all categories in file", [this](std::string_view a)
		{ this->showCategories(a); });
}

MMCQLApplication::~MMCQLApplication()
{
}

void MMCQLApplication::showPagerForData(cif::category &cat)
{
	std::vector<std::string> order;
	for (auto item : cat.get_items())
		order.emplace_back(item);

	std::stringstream os;
	cat.write(os, gOutputFormat, order, false);
	::showPagerForData(os);
}

// --------------------------------------------------------------------
// Some more backslash commands

void MMCQLApplication::loadCifFile(std::string_view f)
{
	if (m_modified and m_file)
		m_file->save(m_filename);

	m_connection.reset();
	m_modified = false;

	try
	{
		m_filename = cif::trim_copy(f);

		cif::gzio::ifstream in(m_filename);
		if (not in.is_open())
			throw std::runtime_error("Could not open file " + m_filename.string());

		m_file.reset(new cif::file{ in });

		if (not m_file->empty())
			loadDatablock(m_file->front().name());
	}
	catch (const std::exception &ex)
	{
		std::cout << "Error loading " << std::quoted(f) << ": " << ex.what() << "\n";
		m_file.reset();
		m_modified = false;
	}
}

void MMCQLApplication::loadDatablock(std::string_view d)
{
	m_dbname = d;

	cif::datablock &db = m_file->operator[](d);
	db.load_dictionary();

	m_connection.reset(new cif::cql::connection(db));
}

void MMCQLApplication::showDatablocks(std::string_view)
{
	if (m_file)
	{
		cif::category datablocks("datablocks");

		for (auto &db : *m_file)
		{
			datablocks.emplace({ //
				{ "name", db.name() } });
		}

		showPagerForData(datablocks);
	}
}

void MMCQLApplication::showCategories(std::string_view)
{
	if (m_file)
	{
		cif::category categories("categories");

		for (auto &db : *m_file)
		{
			for (auto &cat : db)
				categories.emplace({ //
					{ "name", cat.name() },
					{ "datablock", db.name() } });
		}

		showPagerForData(categories);
	}
}

// --------------------------------------------------------------------

void MMCQLApplication::processCommand(std::string_view cmd, std::string_view args)
{
	for (auto &bc : gBackslashCommands)
	{
		if (bc.mCommand != cmd)
			continue;

		try
		{
			bc.mFunc(args);
		}
		catch (const std::exception &ex)
		{
			std::cout << "Error processing command: " << ex.what() << "\n";
		}

		return;
	}

	std::cout << "Unknown command\n";
}

void MMCQLApplication::loop()
{
	using namespace std::literals;
	using namespace readln;

	load_history(".mmcql-history");
	// linenoiseSetMultiLine(1);

	// linenoiseSetCompletionCallback(completion);

	std::cout << "mmql (" << kVersionNumber << ")\n"
			  << "Type \"help\" for help.\n\n";

	std::string line, sql;
	while (getline("cql> ", line))
	{
		add_history(line);

		if (line == "help")
		{
			displayHelp();
			continue;
		}

		// bool executeWithoutSemicolon = false;

		if (line.starts_with('\\'))
		{
			auto args = line.find_first_of(" \t\r\n");
			auto cmd = line.substr(0, args);
			if (cmd == "\\q") // quit
				break;
			// if (cmd == "g")
			// 	executeWithoutSemicolon = true;
			// else
			{
				processCommand(cmd, args == std::string_view::npos ? ""sv : line.substr(args));
				continue;
			}
		}

		if (not m_connection)
		{
			std::cout << "Please load an mmCIF file first using the \\file command\n";
			continue;
		}

		try
		{
			cif::cql::transaction tx(*m_connection);
			auto r = tx.exec(line);

			if (r.empty())
				std::cout << "OK\n";
			else
				showPagerForData(r.get_category());

			tx.commit();
			m_modified = true;
		}
		catch (const std::exception &ex)
		{
			std::cout << "Error executing statement(s): " << ex.what() << "\n";
		}
	}

	if (m_modified and m_file)
		m_file->save(m_filename);

	save_history(".mmcql-history");
}

// -----------------------------------------------------------------------

int pr_main(int argc, char *argv[])
{
	using namespace std::literals;

	auto &config = mcfp::config::instance();

	config.init("mmcql [options] cif-file",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "Verbose output, repeat to increase verbosity level"),

		mcfp::make_option<std::string>("file,f", "Read SQL commands from file"),
		mcfp::make_option<std::string>("command,c", "Single SQL script to execute"),

		mcfp::make_option<std::string>("backup,i", ".bak", "Extension for backup file"),

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

	if (config.has("help") or config.operands().size() != 1)
	{
		std::cerr << config << '\n';
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");

	try
	{
		if (config.has("file") or config.has("command"))
		{
			std::filesystem::path p{ config.operands().front() };
			cif::file file{ p };

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

			std::error_code ec;
			auto backup = p.parent_path() / (p.filename().string() + config.get("backup"));

			if (std::filesystem::exists(backup, ec))
				std::filesystem::remove(backup, ec);

			if (ec)
				std::cerr << "Error removing old backup file: " << ec.message() << '\n';
			else
			{
				std::filesystem::rename(p, backup, ec);
				if (ec)
					std::cerr << "Error creating backup file: " << ec.message() << '\n';
			}

			file.save(p);
		}
		else
		{
			MMCQLApplication app;
			app.loadCifFile(config.operands().front());
			app.loop();
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error in mmcql: " << ex.what() << '\n';
		exit(1);
	}

	return 0;
}
