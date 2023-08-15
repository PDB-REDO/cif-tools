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

#include "revision.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
#include <regex>

#include <mcfp/mcfp.hpp>
#include <cif++.hpp>
#include <cif++/gzio.hpp>

namespace fs = std::filesystem;

class statsParser : public cif::sac_parser
{
  public:
	statsParser(const std::string &file, std::istream &is, const std::string &pattern, bool quiet, bool printLineNr, bool invertMatch)
		: sac_parser(is)
		, mFile(file)
		, mRx(pattern)
		, mQuiet(quiet)
		, mLineNr(printLineNr)
		, mInvertMatch(invertMatch)
	{
	}

	statsParser(const std::string &file, std::istream &is, const std::string &tag, const std::string &pattern, bool quiet, bool printLineNr, bool invertMatch)
		: statsParser(file, is, pattern, quiet, printLineNr, invertMatch)
	{
		std::tie(mCat, mItem) = cif::split_tag_name(tag);
	}

	size_t getMatches() const { return mMatches; }

	void produce_datablock(std::string_view name) override
	{
	}

	void produce_category(std::string_view name) override
	{
	}

	void produce_row() override
	{
	}

	void produce_item(std::string_view category, std::string_view item, std::string_view value) override
	{
		if ((mCat.empty() or cif::iequals(category, mCat)) and
			(mItem.empty() or cif::iequals(item, mItem)) and
			std::regex_search(value.begin(), value.end(), mRx) == not mInvertMatch)
		{
			++mMatches;

			if (not mQuiet)
			{
				if (not mFile.empty())
					std::cout << mFile << ':';
				if (mLineNr)
					std::cout << mLineNr << ':';
				std::cout << value << std::endl;
			}
		}
	}

	std::string mFile;
	std::string mCat, mItem;
	std::regex mRx;
	size_t mMatches = 0;
	bool mQuiet, mLineNr, mInvertMatch;
};

size_t cifGrep(const std::string &pattern, const std::string &tag, const std::string &file, std::istream &is, bool quiet, bool printLineNr, bool invertMatch)
{
	size_t result = 0;

	if (tag.empty())
	{
		statsParser gp(file, is, pattern, quiet, printLineNr, invertMatch);
		gp.parse_file();

		result = gp.getMatches();
	}
	else
	{
		statsParser gp(file, is, tag, pattern, quiet, printLineNr, invertMatch);
		gp.parse_file();

		result = gp.getMatches();
	}

	return result;
}

int pr_main(int argc, char *argv[])
{
	auto &config = mcfp::config::instance();

	config.init(
		"cif-grep [options] pattern file1 [file2...]",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,V", "Verbose output"),

		mcfp::make_option<std::string>("item,i", "The item (tag) to search"),

		mcfp::make_option("quiet,q", "Only print files matching pattern"),
		mcfp::make_option("count,c", "Only show number of hits"),
		mcfp::make_option("invert-match,v", "Select fields NOT matching the pattern"),
		mcfp::make_option("line-number,n", "Print line numbers"),
		mcfp::make_option("no-filename,h", "Don't print the filename"),
		mcfp::make_option("with-filename,H", "Do print the filename"),

		mcfp::make_option("files-with-matches,l", "Print only names of files containing matches"),
		mcfp::make_option("recursive,r", "Search recursively"),

		mcfp::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)")
	);

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().empty())
	{
		std::cerr << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");

	bool quiet = config.has("quiet");
	bool filenamesOnly = config.has("files-with-matches");
	bool countOnly = config.has("count");
	bool noFileNames = filenamesOnly == false and config.has("no-filename");
	bool doFileNames = config.has("with-filename");
	bool lineNumbers = config.has("line-number");
	bool invertMatch = config.has("invert-match");
	size_t count = 0;

	quiet = quiet or countOnly;

	std::string pattern = config.operands().front();
	std::string tag;

	if (config.has("item"))
	{
		tag = config.get<std::string>("item");

		std::string cat, item;
		std::tie(cat, item) = cif::split_tag_name(tag);

		if (cat.empty())
			throw std::runtime_error("Invalid category in tag: '" + cat + '\'');

		if (item.empty())
			throw std::runtime_error("Invalid item: '" + item + '\'');

		if (cif::VERBOSE > 0)
			std::cerr << "matching only for category: " << cat << " and item " << item << std::endl;
	}

	size_t result = false;
	if (config.operands().size() == 1 and not config.has("recursive"))
	{
		result = cifGrep(pattern, tag, "stdin", std::cin, quiet or filenamesOnly, lineNumbers, invertMatch);
		if (doFileNames or (filenamesOnly and result != 0))
			std::cout << "stdin" << std::endl;
		if (countOnly)
			std::cout << result << std::endl;
	}
	else
	{
		std::vector<std::string> files = config.operands();
		files.erase(files.begin());

		if (config.has("recursive"))
		{
			if (files.empty())
				files.push_back(fs::current_path());

			std::vector<std::string> expanded;
			for (auto file : files)
			{
				if (fs::is_directory(file))
				{
					for (auto i = fs::recursive_directory_iterator(file);
						 i != fs::recursive_directory_iterator(); ++i)
					{
						fs::path p = i->path();
						if (fs::is_regular_file(p))
							expanded.push_back(p.string());
					}
				}
				else
					expanded.push_back(file);
			}

			files = expanded;
		}

		std::vector<std::tuple<size_t, std::string>> filesWithSizes;
		size_t totalSize = 0;

		transform(files.begin(), files.end(), back_inserter(filesWithSizes),
			[&totalSize](const std::string &f) -> std::tuple<size_t, std::string>
			{
				size_t size = fs::file_size(f);
				totalSize += size;
				return std::make_tuple(size, f);
			});

		if (doFileNames)
			noFileNames = false;
		else if (files.size() <= 1)
			noFileNames = true;

		for (const auto &[size, file] : filesWithSizes)
		{
			fs::path f(file);

			if (not fs::is_regular_file(f))
				continue;

			if (cif::VERBOSE > 0)
				std::cerr << f << std::endl;

			cif::gzio::ifstream in(f);
			if (not in.is_open())
				throw std::runtime_error("Could not open file " + f.string());

			try
			{
				size_t r = cifGrep(pattern, tag, noFileNames ? "" : f.filename().string(), in, quiet or filenamesOnly, lineNumbers, invertMatch);

				count += r;

				if (cif::VERBOSE or (countOnly and not noFileNames))
					std::cout << f << ':' << r << std::endl;

				if (r > 0)
					result = true;
			}
			catch (const std::exception &e)
			{
				std::cerr << std::endl
						  << "exception for " << f << std::endl
						  << " => " << e.what() << std::endl;
			}
		}
	}

	if (noFileNames and countOnly)
		std::cout << count << std::endl;

	return result ? 0 : 1;
}
