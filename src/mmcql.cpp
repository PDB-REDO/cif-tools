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

#include "cif++/validate.hpp"
#include "mrsrc.hpp"
#include "revision.hpp"

#include <cif++/cif++.hpp>
#include <cif++/category.hpp>
#include <cif++/cql.hpp>
#include <cif++/datablock.hpp>
#include <cif++/gzio.hpp>
#include <cif++/text.hpp>

#include <readln.hpp>

#include <mcfp/mcfp.hpp>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <istream>
#include <memory>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <fcntl.h>
#include <sys/poll.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <system_error>
#include <unistd.h>

#if __has_include(<termios.h>)
#include <termios.h>
#include <sys/ioctl.h>
#endif

// --------------------------------------------------------------------
// Globals to store settings

bool gOutputModeOnlyRows = false;
cif::category::output_format gOutputFormat = cif::category::output_format::column;
std::ofstream gOutputFile;

// --------------------------------------------------------------------
# include <climits>
# include <sys/ioctl.h>
# include <termios.h>

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

std::string to_string(cif::category::output_format fmt)
{
	switch (gOutputFormat)
	{
		case cif::category::output_format::cif:
			return "cif";
		case cif::category::output_format::csv:
			return "csv";
		case cif::category::output_format::tsv:
			return "tsv";
		case cif::category::output_format::list:
			return "list";
		case cif::category::output_format::column:
			return "column";
		case cif::category::output_format::markdown:
			return "markdown";
		case cif::category::output_format::table:
			return "table";
		case cif::category::output_format::box:
			return "box";
		default:
			throw std::runtime_error("Unknown format");
	}
}

void from_string(const std::string &s, cif::category::output_format &fmt)
{
	if (s == "cif")
		fmt = cif::category::output_format::cif;
	else if (s == "csv")
		fmt = cif::category::output_format::csv;
	else if (s == "tsv")
		fmt = cif::category::output_format::tsv;
	else if (s == "list")
		fmt = cif::category::output_format::list;
	else if (s == "column")
		fmt = cif::category::output_format::column;
	else if (s == "markdown")
		fmt = cif::category::output_format::markdown;
	else if (s == "table")
		fmt = cif::category::output_format::table;
	else if (s == "box")
		fmt = cif::category::output_format::box;
	else
		throw std::runtime_error("Unknown output format: " + s);
}

void writeResult(cif::category &cat, std::ostream &os)
{
	cat.write(os, gOutputFormat, {}, true);
}

void writeResult(cif::category &cat)
{
	if (gOutputFile.is_open())
		writeResult(cat, gOutputFile);
	else
		writeResult(cat, std::cout);
}

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
	Formatting,
	DisplayType
};

std::map<CommandCategory, std::string> kCommandCategoryLabels{
	{ CommandCategory::General, "General" },
	{ CommandCategory::Help, "Help" },
	{ CommandCategory::Input_Output, "Input/Output" },
	{ CommandCategory::Formatting, "Formatting" },
	{ CommandCategory::DisplayType, "Display content and type" }
};

struct BackslashCommand
{
	CommandCategory mCategory;
	std::string_view mCommand;
	std::string_view mCommandWithArg;
	std::string_view mDesc;
	std::function<void(std::string_view arg)> mFunc;
};

std::vector<BackslashCommand> gBackslashCommands{
	{ CommandCategory::General,
		"\\copyright", "\\copyright",
		"show copyright and usage", [](std::string_view)
		{
			mrsrc::istream license("LICENSE");
			showPagerForData(license);
		} },
	{ CommandCategory::General, "\\?", "\\?", "show help on backslash commands", [](std::string_view)
		{
			displayOptions();
		} },
	{ CommandCategory::General, "\\h", "\\h", "show help on SQL syntax", [](std::string_view)
		{
			std::cout << "This will eventually be the SQL help page\n";
		} },
	{ CommandCategory::Formatting, "\\format", "\\format [fmt]", "Output format, show current or select one of: cif, csv, tsv, list, column, markdown, table and box", [](std::string_view fmt)
		{
			std::string f{ fmt };
			cif::trim(f);

			if (f.empty())
				std::cout << "Current output format is " << std::quoted(to_string(gOutputFormat)) << "\n";
			else
				from_string(f, gOutputFormat);
		} }
};

void displayOptions()
{
	std::stringstream ss;
	for (const auto &[cat, label] : kCommandCategoryLabels)
	{
		auto cmds = gBackslashCommands | std::views::filter([cat](BackslashCommand &cmd)
											 { return cmd.mCategory == cat; });

		if (cmds.empty())
			continue;

		for (bool first = true; auto cmd : cmds)
		{
			if (std::exchange(first, false))
				ss << label << "\n";
			ss << std::format("  {:15}  {}\n", cmd.mCommandWithArg, cmd.mDesc);
		}
		ss << "\n";
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

	~MMCQLApplication();

  private:
	void loadCifFile(std::string_view f);
	void loadDictionary(std::string_view d);
	void loadDatablock(std::string_view d);

	void showPagerForData(cif::category &cat);

	void showCategories(std::string_view, bool show_all);
	void showDatablocks(std::string_view);
	void showItems(std::string_view, bool with_aliases);

	void processCommand(std::string_view cmd, std::string_view args);

	std::filesystem::path m_file_name;
	std::unique_ptr<cif::file> m_file;
	std::string m_db_name, m_dict_name;
	std::unique_ptr<cif::cql::connection> m_connection;
};

// --------------------------------------------------------------------

MMCQLApplication::MMCQLApplication()
{
	auto &config = mcfp::config::instance();

	if (config.operands().size() > 0)
		m_file_name = config.operands().front();
	if (config.has("datablock"))
		m_db_name = config.get("datablock");
	if (config.has("dict"))
		m_dict_name = config.get("dict");

	gBackslashCommands.emplace_back(
		CommandCategory::Input_Output,
		"\\file", "\\file <file>", "load new mmCIF file", [this](std::string_view f)
		{ this->loadCifFile(f); });
	gBackslashCommands.emplace_back(
		CommandCategory::Input_Output,
		"\\db", "\\db <datablock>", "load a datablock from the current file", [this](std::string_view f)
		{ this->loadDatablock(f); });
	gBackslashCommands.emplace_back(
		CommandCategory::Input_Output,
		"\\dict", "\\dict <dictionary>", "use a specified dictionary for the current datablock", [this](std::string_view f)
		{ this->loadDictionary(f); });
	gBackslashCommands.emplace_back(
		CommandCategory::DisplayType,
		"\\dd", "\\dd", "show datablocks in file", [this](std::string_view a)
		{ this->showDatablocks(a); });
	gBackslashCommands.emplace_back(
		CommandCategory::DisplayType,
		"\\dc", "\\dc", "show all categories in current datablock", [this](std::string_view a)
		{ this->showCategories(a, false); });
	gBackslashCommands.emplace_back(
		CommandCategory::DisplayType,
		"\\dc+", "\\dc+", "show all categories in current file", [this](std::string_view a)
		{ this->showCategories(a, true); });
	gBackslashCommands.emplace_back(
		CommandCategory::DisplayType,
		"\\d", "\\d <category>", "show items in category", [this](std::string_view a)
		{ this->showItems(a, false); });
	gBackslashCommands.emplace_back(
		CommandCategory::DisplayType,
		"\\d+", "\\d+ <category>", "show items in category including aliased names", [this](std::string_view a)
		{ this->showItems(a, true); });
}

MMCQLApplication::~MMCQLApplication()
{
	if (m_connection and m_connection->is_modified() and m_file)
		m_file->save(m_file_name);
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
	if (m_connection and m_connection->is_modified() and m_file)
		m_file->save(m_file_name);

	m_connection.reset();

	try
	{
		m_file_name = cif::trim_copy(f);

		cif::gzio::ifstream in(m_file_name);
		if (not in.is_open())
			throw std::runtime_error("Could not open file " + m_file_name.string());

		m_file.reset(new cif::file{ in });

		if (not m_file->empty())
		{
			if (m_file->contains(m_db_name))
				loadDatablock(m_db_name);
			else
				loadDatablock(m_file->front().name());
		}

		std::cout << "Loaded file " << m_file_name << " and datablock " << m_db_name << "\n";
	}
	catch (const std::exception &ex)
	{
		std::cout << "Error loading " << std::quoted(f) << ": " << ex.what() << "\n";
		m_file.reset();
	}
}

void MMCQLApplication::loadDatablock(std::string_view d)
{
	m_db_name = cif::trim_copy(d);

	cif::datablock &db = m_file->operator[](m_db_name);

	if (m_dict_name.empty())
		db.load_dictionary();
	else
		db.load_dictionary(m_dict_name);

	m_connection.reset(new cif::cql::connection(db));
}

void MMCQLApplication::loadDictionary(std::string_view dn)
{
	m_dict_name = cif::trim_copy(dn);

	cif::datablock &db = m_file->operator[](m_db_name);

	if (m_dict_name.empty())
	{
		cif::category audit_conform;

		if (auto v = db.get_validator(); v != nullptr)
			v->fill_audit_conform(audit_conform);

		if (audit_conform.empty())
			std::cout << "No dictionary loaded\n";
		else
			showPagerForData(audit_conform);
	}
	else
	{
		db.load_dictionary(m_dict_name);
		m_connection.reset(new cif::cql::connection(db));
	}
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

void MMCQLApplication::showCategories(std::string_view, bool show_all)
{
	if (m_file)
	{
		cif::category categories("categories");

		if (show_all)
		{
			for (auto &db : *m_file)
			{
				for (auto &cat : db)
					categories.emplace({ //
						{ "name", cat.name() },
						{ "datablock", db.name() } });
			}
		}
		else
		{
			for (auto &db : *m_file)
			{
				if (db.name() != m_db_name)
					continue;
				for (auto &cat : db)
					categories.emplace({ //
						{ "name", cat.name() } });
			}
		}

		showPagerForData(categories);
	}
}

void MMCQLApplication::showItems(std::string_view a, bool with_aliases)
{
	std::string name = cif::trim_copy(a);

	if (m_file and not m_db_name.empty())
	{
		cif::category items("items in " + std::string{ name });

		for (auto &cat : m_file->operator[](m_db_name))
		{
			if (cat.name() != name)
				continue;

			if (auto cv = cat.get_cat_validator(); cv != nullptr)
			{
				for (auto iv : cv->m_item_validators)
				{
					auto ri = items.emplace({ { "name", iv.m_item_name },
						{ "type", iv.m_type->m_name },
						{ "mandatory", iv.m_mandatory },
						{ "default", iv.m_default },
						{ "enumeration", cif::join(iv.m_enums, ",") } });

					if (with_aliases and not iv.m_aliases.empty())
					{
						std::ostringstream os;
						for (bool first = true; auto &alias : iv.m_aliases)
						{
							if (std::exchange(first, false) == false)
								os << ";";
							os << alias.m_name << " in " << alias.m_dict << '@' << alias.m_vers;
						}

						ri->assign("aliases", os.str(), false, false);
					}
				}
			}
			else
			{
				for (auto item : cat.get_items())
					items.emplace({ { "name", item } });
			}
		}

		showPagerForData(items);
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

	if (not m_file_name.empty())
		loadCifFile(m_file_name.string());

	std::string line, sql;
	while (getline("cql> ", line))
	{
		cif::trim(line);
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

			processCommand(cmd, args == std::string_view::npos ? ""sv : line.substr(args));
			continue;
		}

		if (not m_connection)
		{
			std::cout << "Please load an mmCIF file first using the \\file command\n";
			continue;
		}

		if (sql.empty())
			sql = line;
		else if (std::isspace(sql.back()))
			sql += line;
		else
			sql = sql + ' ' + line;

		while (m_connection->is_complete_statement(sql))
		{
			try
			{
				auto r = m_connection->exec(sql, sql);

				if (r.empty())
					std::cout << "OK\n";
				else
					showPagerForData(r.get_category());
			}
			catch (const std::exception &ex)
			{
				sql.clear();
				std::cout << "Error executing statement(s): " << ex.what() << "\n";
			}
		}
	}

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

		mcfp::make_option<std::string>("datablock,D", "Datablock to use, default is first"),
		mcfp::make_option<std::string>("dict,d", "Use specified mmcif dictionary"),

		mcfp::make_option<std::string>("file,f", "Read SQL commands from file"),
		mcfp::make_option<std::string>("command,c", "Single SQL script to execute"),

		mcfp::make_option<std::string>("format", "column", "Output format to use"));

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

	from_string(config.get("format"), gOutputFormat);

	try
	{
		if (config.has("file") or config.has("command"))
		{
			std::filesystem::path p{ config.operands().front() };
			cif::file file{ p };

			cif::datablock &db = config.has("data-block") ? file[config.get("data-block")] : file.front();

			if (config.has("dict"))
				db.load_dictionary(config.get("dict"));
			else
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
					writeResult(r.get_category());
				tx.commit();
			}
			else if (config.has("command"))
			{
				cif::cql::transaction tx(connection);
				auto r = tx.exec(config.get("command"));
				if (not r.empty())
					writeResult(r.get_category());
				tx.commit();
			}

			if (connection.is_modified())
			{
				std::error_code ec;
				auto backup = p.parent_path() / (p.filename().string() + ".bak");

				if (std::filesystem::exists(backup, ec))
					std::filesystem::remove(backup, ec);

				if (ec)
					std::cerr << "Error removing old backup file: " << ec.message() << '\n';

				std::filesystem::rename(p, backup, ec);
				if (ec)
					std::cerr << "Error creating backup file: " << ec.message() << '\n';

				file.save(p);
			}
		}
		else
		{
			MMCQLApplication app;
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
