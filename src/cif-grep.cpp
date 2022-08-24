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

#include <filesystem>
#include <fstream>
#include <functional>
#include <regex>

#include <gxrio.hpp>

#include <cif++.hpp>

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

	void produce_datablock(const std::string &name) override
	{
	}

	void produce_category(const std::string &name) override
	{
	}

	void produce_row() override
	{
	}

	void produce_item(const std::string &category, const std::string &item, const std::string &value) override
	{
		if ((mCat.empty() or cif::iequals(category, mCat)) and
			(mItem.empty() or cif::iequals(item, mItem)) and
			std::regex_search(value, mRx) == not mInvertMatch)
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
	po::options_description visible_options("cif-grep [option...] pattern [file ...]");
	visible_options.add_options()("item,i", po::value<std::string>(), "Item tag to scan, default is all item values")("help", "Display help message")("version", "Print version")("quiet,q", "Only print files matching pattern")("count,c", "Only show number of hits")("invert-match,v", "Select fields NOT matching the pattern")("line-number,n", "Print line numbers")("no-filename,h", "Don't print the filename")("with-filename,H", "Do print the filename")("verbose,V", "Verbose output")("files-with-matches,l", "Print only names of files containing matches")("recursive,r", "Search recursively");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()("pattern", po::value<std::string>(), "Pattern")("input", po::value<std::vector<std::string>>(), "Input files")("debug,d", po::value<int>(), "Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("pattern", 1);
	p.add("input", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		write_version_string(std::cout, vm.count("verbose"));
		exit(0);
	}

	if (vm.count("help") or vm.count("pattern") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(vm.count("help") ? 0 : 1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	bool quiet = vm.count("quiet") > 0;
	bool filenamesOnly = vm.count("files-with-matches") > 0;
	bool countOnly = vm.count("count") > 0;
	bool noFileNames = filenamesOnly == false and vm.count("no-filename") > 0;
	bool doFileNames = vm.count("with-filename") > 0;
	bool lineNumbers = vm.count("line-number") > 0;
	bool invertMatch = vm.count("invert-match") > 0;
	size_t count = 0;

	quiet = quiet or countOnly;

	std::string pattern = vm["pattern"].as<std::string>();
	std::string tag;
	if (vm.count("item"))
	{
		tag = vm["item"].as<std::string>();

		std::string cat, item;
		std::tie(cat, item) = cif::split_tag_name(tag);

		if (cat.empty())
			throw std::runtime_error("Invalid category in tag: '" + cat + '\'');

		if (item.empty())
			throw std::runtime_error("Invalid item: '" + item + '\'');

		if (cif::VERBOSE)
			std::cerr << "matching only for category: " << cat << " and item " << item << std::endl;
	}

	size_t result = false;
	if (vm.count("input") == 0 and not vm.count("recursive"))
	{
		result = cifGrep(pattern, tag, "stdin", std::cin, quiet or filenamesOnly, lineNumbers, invertMatch);
		if (doFileNames or (filenamesOnly and result != 0))
			std::cout << "stdin" << std::endl;
		if (countOnly)
			std::cout << result << std::endl;
	}
	else
	{
		std::vector<std::string> files;
		if (vm.count("input"))
			files = vm["input"].as<std::vector<std::string>>();

		if (vm.count("recursive"))
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

			if (cif::VERBOSE)
				std::cerr << f << std::endl;

			gxrio::ifstream in(f);
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
