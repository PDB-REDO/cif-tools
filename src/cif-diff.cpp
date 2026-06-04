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

#include "cif++/category.hpp"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <ranges>
#include <stdexcept>
#include <system_error>
#ifndef WIN32
# include <sys/wait.h>
#endif

#include <filesystem>
#include <fstream>
#include <functional>

#include <cif++/cif++.hpp>
#include <mcfp/mcfp.hpp>
#include <utility>

#include "revision.hpp"

using cif::row_handle;

namespace fs = std::filesystem;

class fd_streambuf : public std::streambuf
{
  public:
	using int_type = std::streambuf::int_type;
	using char_type = std::streambuf::char_type;
	using traits_type = std::streambuf::traits_type;

	fd_streambuf(int fd)
		: m_fd(fd)
	{
	}

	~fd_streambuf() override
	{
		close();
	}

	void close()
	{
		if (m_fd >= 0)
		{
			sync();
			::close(m_fd);
		}

		m_fd = -1;
	}

	int_type overflow(int_type ch) override
	{
		assert(pptr() == epptr());

		if (m_fd < 0)
			return traits_type::eof();

		if (ch == traits_type::eof() or sync() == -1)
			return traits_type::eof();

		*pptr() = traits_type::to_char_type(ch);
		pbump(1);
		return ch;
	}

	int sync() override
	{
		int result = 0;

		if (m_fd >= 0)
		{
			char *p = pbase();
			while (p < pptr() and result == 0)
			{
				int n = ::write(m_fd, p, pptr() - p);
				if (n <= 0)
					result = -1;
				p += n;
			}

			if (result == 0)
				setp(m_buffer.data(), m_buffer.data() + m_buffer.size());
		}

		return result;
	}

  private:
	static const size_t kBufferSize = 4096;

	int m_fd;
	std::array<char_type, kBufferSize> m_buffer;
};

// --------------------------------------------------------------------

class templateParser : public cif::sac_parser
{
  public:
	templateParser(std::istream &is)
		: sac_parser(is)
	{
	}

	void produce_datablock(std::string_view name) override
	{
	}

	void produce_category(std::string_view name) override
	{
	}

	void produce_row() override
	{
	}

	void produce_item(std::string_view category, std::string_view item, cif::item_value value) override
	{
		std::ostringstream tag;
		tag << '_' << category << '.' << item;
		if (std::ranges::find(std::views::reverse(mOrder), tag.str()) == mOrder.rend())
			mOrder.push_back(tag.str());
	}

	std::vector<std::string> mOrder;
};

// --------------------------------------------------------------------

void compareCategories(cif::category &a, cif::category &b, size_t maxDiffCount)
{
	using namespace std::placeholders;

	auto validator = a.get_validator();
	auto catValidator = validator->get_validator_for_category(a.name());
	if (catValidator == nullptr)
		throw std::runtime_error("missing cat validator");

	using compType = std::function<int(std::string_view, std::string_view)>;
	std::vector<std::tuple<std::string, compType>> tags;
	auto keys = catValidator->m_keys;

	for (auto &item : a.key_items())
	{
		auto iv = catValidator->get_validator_for_item(item);
		if (iv == nullptr)
			throw std::runtime_error("missing item validator");
		auto tv = iv->m_type;
		if (tv == nullptr)
			throw std::runtime_error("missing type validator");
		tags.emplace_back(item,
			[tv](std::string_view a, std::string_view b)
			{ return tv->compare(a, b); });

		auto pred = [item](const std::string &s) -> bool
		{
			return cif::iequals(item, s) == 0;
		};
	}

	a.reorder_by_index();
	b.reorder_by_index();

	auto rowLess = [&](const cif::row_handle &a, const cif::row_handle &b) -> bool
	{
		int d = 0;

		for (const auto &[tag, compare] : tags)
		{
			d = a[tag].compare(b[tag]);

			if (d != 0)
				break;
		}

		return d < 0;
	};

	//	std::vector<cif::row_handle> rowsA(a.begin(), a.end()), rowsB(b.begin(), b.end());
	//	sort(rowsA.begin(), rowsA.end(), rowLess);
	//	sort(rowsB.begin(), rowsB.end(), rowLess);

	auto ai = a.begin(), bi = b.begin();

	struct Diff
	{
		virtual ~Diff() = default;

		std::string key(cif::row_handle r, std::vector<std::string> &keys)
		{
			std::vector<std::string> v;
			for (auto k : keys)
				v.push_back(r[k].get<std::string>());
			return "[" + cif::join(v, ", ") + "]";
		}

		virtual void report(std::vector<std::string> &keys) = 0;
	};

	struct ExtraADiff : public Diff
	{
		cif::row_handle A;

		ExtraADiff(cif::row_handle r)
			: A(std::move(r))
		{
		}

		void report(std::vector<std::string> &keys) override
		{
			std::cout << "Extra row in A with key " << key(A, keys) << '\n';
		}
	};

	struct ExtraBDiff : public Diff
	{
		cif::row_handle B;

		ExtraBDiff(cif::row_handle r)
			: B(std::move(r))
		{
		}

		void report(std::vector<std::string> &keys) override
		{
			std::cout << "Extra row in B with key " << key(B, keys) << '\n';
		}
	};

	struct ValueDiff : public Diff
	{
		cif::row_handle A, B;
		std::vector<std::string> missingA, missingB, different;

		ValueDiff(cif::row_handle a, row_handle b, std::vector<std::string> &&missingA, std::vector<std::string> &&missingB, std::vector<std::string> &&different)
			: A(std::move(a))
			, B(std::move(b))
			, missingA(std::move(missingA))
			, missingB(std::move(missingB))
			, different(std::move(different))
		{
		}

		void report(std::vector<std::string> &keys) override
		{
			std::cout << "Differences in rows with key " << key(A, keys) << '\n';

			for (auto &item : different)
			{
				std::cout << "    " << item << " (A): '" << A[item].get<std::string>() << '\'' << '\n'
						  << "    " << item << " (B): '" << B[item].get<std::string>() << '\'' << '\n';
			}

			for (auto &item : missingA)
			{
				std::cout << "    " << item << " (A): <missing>\n"
						  << "    " << item << " (B): '" << B[item].get<std::string>() << '\'' << '\n';
			}

			for (auto &item : missingB)
			{
				std::cout << "    " << item << " (A): '" << A[item].get<std::string>() << '\'' << '\n'
						  << "    " << item << " (B): <missing>\n";
			}
		}
	};

	std::vector<Diff *> diffs;

	while ((maxDiffCount == 0 or diffs.size() < maxDiffCount) and (ai != a.end() or bi != b.end()))
	{
		if (ai == a.end())
		{
			diffs.push_back(new ExtraBDiff{ *bi++ });
			continue;
		}

		if (bi == b.end())
		{
			diffs.push_back(new ExtraADiff{ *ai++ });
			continue;
		}

		cif::row_handle ra = *ai, rb = *bi;

		if (rowLess(ra, rb))
		{
			diffs.push_back(new ExtraADiff{ *ai++ });
			continue;
		}

		if (rowLess(rb, ra))
		{
			diffs.push_back(new ExtraBDiff{ *bi++ });
			continue;
		}

		std::vector<std::string> missingA, missingB, different;

		for (auto &tt : tags)
		{
			std::string tag;
			compType compare;

			tie(tag, compare) = tt;

			// make it an option to compare unapplicable to empty or something
			if (ra[tag].empty())
				missingA.push_back(tag);
			else if (rb[tag].empty())
				missingB.push_back(tag);
			else if (ra[tag].compare(rb[tag]) != 0)
				different.push_back(tag);
		}

		++ai;
		++bi;

		if (not missingA.empty() or not missingB.empty() or not different.empty())
			diffs.push_back(new ValueDiff{ ra, rb, std::move(missingA), std::move(missingB), std::move(different) });
	}

	if (not diffs.empty())
	{
		std::cout << std::string(cif::get_terminal_width(), '-') << '\n'
				  << "Differences in values for category " << a.name() << '\n'
				  << '\n';

		for (auto diff : diffs)
		{
			diff->report(keys);
			delete diff;
		}

		if (diffs.size() == maxDiffCount)
			std::cout << "...\n";
		std::cout << '\n';
	}
}

void compareCifs(cif::datablock &dbA, cif::datablock &dbB, const cif::iset &categories, int maxDiffCount)
{
	if (dbA.empty() or dbB.empty())
		throw std::runtime_error("Invalid emtpy datablocks");

	std::vector<std::string> catA, catB;

	for (auto &cat : dbA)
		catA.push_back(cat.name());
	std::ranges::sort(catA, [](const std::string &a, const std::string &b)
		{ return cif::icompare(a, b) < 0; });

	for (auto &cat : dbB)
		catB.push_back(cat.name());
	std::ranges::sort(catB, [](const std::string &a, const std::string &b)
		{ return cif::icompare(a, b) < 0; });

	// loop over categories twice, to group output
	// First iteration is to list missing categories.

	std::vector<std::string> missingA, missingB;

	auto catA_i = catA.begin(), catB_i = catB.begin();

	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		cif::to_lower(nA);

		std::string nB = *catB_i;
		cif::to_lower(nB);

		int d = nA.compare(nB);
		if (d > 0)
		{
			auto &cat = dbB[*catB_i++];

			if (not cat.empty())
				missingA.push_back(cat.name());
		}
		else if (d < 0)
		{
			auto &cat = dbA[*catA_i++];

			if (not cat.empty())
				missingB.push_back(cat.name());
		}
		else
			++catA_i, ++catB_i;
	}

	while (catA_i != catA.end())
		missingB.push_back(*catA_i++);

	while (catB_i != catB.end())
		missingA.push_back(*catB_i++);

	if (categories.empty())
	{
		if (not missingA.empty())
			std::cout << "Categories missing in A: " << cif::join(missingA, ", ") << '\n'
					  << '\n';

		if (not missingB.empty())
			std::cout << "Categories missing in B: " << cif::join(missingB, ", ") << '\n'
					  << '\n';
	}

	// Second loop, now compare category values
	catA_i = catA.begin(), catB_i = catB.begin();

	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		cif::to_lower(nA);

		std::string nB = *catB_i;
		cif::to_lower(nB);

		int d = nA.compare(nB);
		if (d > 0)
			++catB_i;
		else if (d < 0)
			++catA_i;
		else
		{
			if (categories.empty() or categories.count(nA))
				compareCategories(dbA[*catA_i], dbB[*catB_i], maxDiffCount);
			++catA_i;
			++catB_i;
		}
	}
}

#ifndef WIN32
void compareCifsText(const std::string &editor, cif::file &a, cif::file &b, std::filesystem::path file_a, std::filesystem::path file_b, bool icase, bool iwhite)
{
	// temp files for external diff

	std::string dir_s = (fs::temp_directory_path() / "cif-diff-XXXXXX").string();
	if (mkdtemp(dir_s.data()) == nullptr)
		throw std::system_error(std::error_code(errno, std::system_category()), "Error creating temporary directory");

	std::filesystem::path dir(dir_s);

	auto out_1 = dir / file_a.filename();
	auto out_2 = dir / file_b.filename();

	if (out_1.extension() == ".gz")
		out_1.replace_extension();
	if (out_2.extension() == ".gz")
		out_2.replace_extension();

	if (out_1 == out_2)
	{
		out_1.replace_filename(out_1.filename().replace_extension(".1" + out_1.extension().string()));
		out_2.replace_filename(out_2.filename().replace_extension(".2" + out_2.extension().string()));
	}

	std::ofstream f1(out_1);
	std::ofstream f2(out_2);

	if (not(f1.is_open() and f2.is_open()))
		throw std::runtime_error("Could not open files for output");

	auto dia = a.begin();
	auto dib = b.begin();

	while (dia != a.end() and dib != b.end())
	{
		auto &da = *dia++;
		auto &db = *dib++;

		f1 << "data_" << da.name() << "\n# \n";
		f2 << "data_" << db.name() << "\n# \n";

		std::vector<std::string> categories;
		for (auto &cat : da)
			categories.emplace_back(cat.name());
		for (auto &cat : db)
		{
			if (not std::ranges::contains(categories, cat.name()))
				categories.emplace_back(cat.name());
		}

		for (auto &category : categories)
		{
			auto cat_a = da.get(category);
			auto cat_b = db.get(category);

			if (not cat_a)
			{
				cat_b->write(f1);
				continue;
			}
			
			if (not cat_b)
			{
				cat_a->write(f1);
				continue;
			}

			cat_a->drop_empty_items();
			cat_b->drop_empty_items();

			std::vector<std::string> items = cat_a->get_items();

			for (auto &item : cat_b->get_items())
			{
				if (not std::ranges::contains(items, item))
					items.emplace_back(item);
			}

			cat_a->write(f1, items, true);
			cat_b->write(f2, items, true);
		}
	}

	while (dia != a.end())
	{
		auto &da = *dia++;
		da.write(f1);
	}

	while (dib != b.end())
	{
		auto &db = *dib++;
		db.write(f2);
	}

	f1.close();
	f2.close();

	std::ostringstream cmd;
	cmd << editor << " -d";

	if (icase)
		cmd << " -c 'set diffopt+=icase'";

	if (iwhite)
		cmd << " -c 'set diffopt-=iwhite'";

	cmd << " " << std::quoted(out_1.string()) << " " << std::quoted(out_2.string());

	switch (auto pid = fork())
	{
		case -1:
			std::cerr << "fork failed: " << std::error_code(errno, std::system_category()).message() << "\n";
			break;
		case 0:
			execlp("/bin/sh", "sh", "-c", cmd.str().c_str(), nullptr);
			std::cerr << "exec of editor failed: " << std::error_code(errno, std::system_category()).message() << "\n";
			exit(-1);
			break;
		default:
			waitpid(pid, nullptr, 0);
			break;
	}

	// Only with vimdiff it is safe to remove the files now
	if (editor == "vim" or editor == "vimdiff")
		std::filesystem::remove_all(dir);
}
#endif

int pr_main(int argc, char *argv[])
{
	auto &config = mcfp::config::instance();

	config.init(
		"cif-diff [options] file1 file2",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "Verbose output"),
		mcfp::make_option<std::vector<std::string>>("category", "Limit comparison to this category, default is all categories. Can be specified multiple times"),
		mcfp::make_option<int>("max-diff-count", 5, "Maximum number of diff items per category, enter zero (0) for unlimited, default is 5"),
		mcfp::make_option<std::string>("editor", "vim", "Editor to use for showing the textual differences. Default is vim, alternative is 'terminal' to dump to stdout."),
		mcfp::make_option("icase", "Ignore case (vimdiff option)"),
		mcfp::make_option("iwhite", "Ignore whitespace (vimdiff option)"),
		mcfp::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)"));

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().size() != 2)
	{
		std::cerr << config << '\n';
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");

	int maxDiffCount = config.get<int>("max-diff-count");

	cif::iset categories;
	if (config.has("category"))
	{
		for (auto cs : config.get<std::vector<std::string>>("category"))
		{
			for (auto cat : cif::split(cs, ",; ", true))
				categories.emplace(cif::to_lower_copy(cat));
		}
	}

	auto input = config.operands();

	cif::gzio::ifstream if1{ input[0] };
	if (not if1.is_open())
		throw std::runtime_error("Could not open file " + input[0]);

	cif::gzio::ifstream if2(input[1]);
	if (not if2.is_open())
		throw std::runtime_error("Could not open file " + input[1]);

	cif::file file1{ if1 };
	cif::file file2{ if2 };

	if (config.get("editor") == "terminal")
		compareCifs(file1.front(), file2.front(), categories, maxDiffCount);
	else
		compareCifsText(config.get("editor"), file1, file2, fs::path(input[0]), fs::path(input[1]), config.has("icase"), config.has("iwhite"));

	return 0;
}
