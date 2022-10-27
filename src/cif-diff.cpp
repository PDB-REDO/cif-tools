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

#include <sys/wait.h>

#include <filesystem>
#include <fstream>
#include <functional>

#include <cif++.hpp>
#include <cfg.hpp>
#include <gxrio.hpp>

#include "revision.hpp"

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

	~fd_streambuf()
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

	virtual int_type overflow(int_type ch)
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

	int sync()
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
				setp(m_buffer.begin(), m_buffer.end());
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
		std::string tag = "_" + category + "." + item;
		if (find(mOrder.rbegin(), mOrder.rend(), tag) == mOrder.rend())
			mOrder.push_back(tag);
	}

	std::vector<std::string> mOrder;
};

// --------------------------------------------------------------------

void compareCategories(cif::category &a, cif::category &b, size_t maxDiffCount)
{
	using namespace std::placeholders;

	//	set<std::string> tagsA(a.fields()), tagsB(b.fields());
	//
	//	if (tagsA != tagsB)
	//		std::cout << "Unequal number of fields" << std::endl;

	auto validator = a.get_validator();
	auto catValidator = validator->get_validator_for_category(a.name());
	if (catValidator == nullptr)
		throw std::runtime_error("missing cat validator");

	typedef std::function<int(std::string_view, std::string_view)> compType;
	std::vector<std::tuple<std::string, compType>> tags;
	auto keys = catValidator->m_keys;
	std::vector<size_t> keyIx;

	for (auto &tag : a.fields())
	{
		auto iv = catValidator->get_validator_for_item(tag);
		if (iv == nullptr)
			throw std::runtime_error("missing item validator");
		auto tv = iv->m_type;
		if (tv == nullptr)
			throw std::runtime_error("missing type validator");
		tags.push_back(std::make_tuple(tag, std::bind(&cif::type_validator::compare, tv, std::placeholders::_1, std::placeholders::_2)));

		auto pred = [tag](const std::string &s) -> bool
		{ return cif::iequals(tag, s) == 0; };
		if (find_if(keys.begin(), keys.end(), pred) == keys.end())
			keyIx.push_back(tags.size() - 1);
	}

	a.reorder_by_index();
	b.reorder_by_index();

	auto rowLess = [&](const cif::row_handle &a, const cif::row_handle &b) -> bool
	{
		int d = 0;

		for (auto kix : keyIx)
		{
			std::string tag;
			compType compare;

			tie(tag, compare) = tags[kix];

			d = compare(a[tag].text(), b[tag].text());

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
		virtual ~Diff() {}

		std::string key(cif::row_handle r, std::vector<std::string> &keys)
		{
			std::vector<std::string> v;
			for (auto k : keys)
				v.push_back(r[k].as<std::string>());
			return "[" + cif::join(v, ", ") + "]";
		}

		virtual void report(std::vector<std::string> &keys) = 0;
	};

	struct ExtraADiff : public Diff
	{
		cif::row_handle A;

		ExtraADiff(cif::row_handle r)
			: A(r)
		{
		}

		virtual void report(std::vector<std::string> &keys)
		{
			std::cout << "Extra row in A with key " << key(A, keys) << std::endl;
		}
	};

	struct ExtraBDiff : public Diff
	{
		cif::row_handle B;

		ExtraBDiff(cif::row_handle r)
			: B(r)
		{
		}

		virtual void report(std::vector<std::string> &keys)
		{
			std::cout << "Extra row in B with key " << key(B, keys) << std::endl;
		}
	};

	struct ValueDiff : public Diff
	{
		cif::row_handle A, B;
		std::vector<std::string> missingA, missingB, different;

		ValueDiff(cif::row_handle a, cif::row_handle b, std::vector<std::string> &&missingA, std::vector<std::string> &&missingB, std::vector<std::string> &&different)
			: A(a)
			, B(b)
			, missingA(move(missingA))
			, missingB(move(missingB))
			, different(move(different))
		{
		}

		virtual void report(std::vector<std::string> &keys)
		{
			std::cout << "Differences in rows with key " << key(A, keys) << std::endl;

			for (auto &item : different)
			{
				std::cout << "    " << item << " (A): '" << A[item].as<std::string>() << '\'' << std::endl
						  << "    " << item << " (B): '" << B[item].as<std::string>() << '\'' << std::endl;
			}

			for (auto &item : missingA)
			{
				std::cout << "    " << item << " (A): <missing>" << std::endl
						  << "    " << item << " (B): '" << B[item].as<std::string>() << '\'' << std::endl;
			}

			for (auto &item : missingB)
			{
				std::cout << "    " << item << " (A): '" << A[item].as<std::string>() << '\'' << std::endl
						  << "    " << item << " (B): <missing>" << std::endl;
			}
		}
	};

	std::vector<Diff *> diffs;

	while ((maxDiffCount == 0 or diffs.size() < maxDiffCount) and (ai != a.end() or bi != b.end()))
	{
		if (ai == a.end())
		{
			diffs.push_back(new ExtraBDiff{*bi++});
			continue;
		}

		if (bi == b.end())
		{
			diffs.push_back(new ExtraADiff{*ai++});
			continue;
		}

		cif::row_handle ra = *ai, rb = *bi;

		if (rowLess(ra, rb))
		{
			diffs.push_back(new ExtraADiff{*ai++});
			continue;
		}

		if (rowLess(rb, ra))
		{
			diffs.push_back(new ExtraBDiff{*bi++});
			continue;
		}

		std::vector<std::string> missingA, missingB, different;

		for (auto &tt : tags)
		{
			std::string tag;
			compType compare;

			tie(tag, compare) = tt;

			// make it an option to compare unapplicable to empty or something

			std::string_view ta = ra[tag].text();
			if (ta == ".")
				ta = "";
			std::string_view tb = rb[tag].text();
			if (tb == ".")
				tb = "";

			if (compare(ta, tb) != 0)
			{
				if (ta.empty())
					missingA.push_back(tag);
				else if (tb.empty())
					missingB.push_back(tag);
				else
					different.push_back(tag);
			}
		}

		++ai;
		++bi;

		if (not missingA.empty() or not missingB.empty() or not different.empty())
			diffs.push_back(new ValueDiff{ra, rb, move(missingA), move(missingB), move(different)});
	}

	if (not diffs.empty())
	{
		std::cout << std::string(cif::get_terminal_width(), '-') << std::endl
				  << "Differences in values for category " << a.name() << std::endl
				  << std::endl;

		for (auto diff : diffs)
		{
			diff->report(keys);
			delete diff;
		}

		if (diffs.size() == maxDiffCount)
			std::cout << "..." << std::endl;
		std::cout << std::endl;
	}
}

void compareCifs(cif::datablock &dbA, cif::datablock &dbB, const cif::iset &categories, int maxDiffCount)
{
	if (dbA.empty() or dbB.empty())
		throw std::runtime_error("Invalid emtpy datablocks");

	std::vector<std::string> catA, catB;

	for (auto &cat : dbA)
		catA.push_back(cat.name());
	sort(catA.begin(), catA.end());

	for (auto &cat : dbB)
		catB.push_back(cat.name());
	sort(catB.begin(), catB.end());

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
			std::cout << "Categories missing in A: " << cif::join(missingA, ", ") << std::endl
					  << std::endl;

		if (not missingB.empty())
			std::cout << "Categories missing in B: " << cif::join(missingB, ", ") << std::endl
					  << std::endl;
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

void compareCifsText(cif::file &a, cif::file &b, const std::string &name_a, const std::string &name_b, bool icase, bool iwhite)
{
	// temp files for vimdiff

	fs::path file_a(name_a);
	fs::path file_b(name_b);

	std::string generated = fs::temp_directory_path() / ("cif-diff-" + file_a.filename().string() + "-XXXXXX.cif");
	std::string original = fs::temp_directory_path() / ("cif-diff-" + file_b.filename().string() + "-XXXXXX.cif");

	int fd[2];

	if ((fd[0] = mkstemps(generated.data(), 4)) < 0 or (fd[1] = mkstemps(original.data(), 4)) < 0)
	{
		std::cerr << "Error creating temp files:  " << strerror(errno) << std::endl;
		exit(1);
	}

	{
		fd_streambuf sb(fd[0]);
		std::ostream out(&sb);
		a.front().write(out);
	}

	// Next the converted cif file

	{
		fd_streambuf sb(fd[1]);
		std::ostream out(&sb);

		b.front().write(out, a.front().get_tag_order());
	}

	std::vector<const char *> nArgv = {
		"/usr/bin/vimdiff"};

	if (icase)
	{
		nArgv.push_back("-c");
		nArgv.push_back("set diffopt+=icase");
	}

	if (iwhite)
	{
		nArgv.push_back("-c");
		nArgv.push_back("set diffopt-=iwhite");
	}

	nArgv.push_back(original.c_str());
	nArgv.push_back(generated.c_str());
	nArgv.push_back(nullptr);

	int pid = fork();

	if (pid <= 0)
	{
		if (execv(nArgv[0], const_cast<char *const *>(nArgv.data())) < 0)
			std::cerr << "Failed to execute vimdiff" << std::endl;
		exit(1);
	}

	int status;
	waitpid(pid, &status, 0);

	if (WIFEXITED(status))
	{
		unlink(generated.c_str());
		unlink(original.c_str());
	}
}

int pr_main(int argc, char *argv[])
{
	auto &config = cfg::config::instance();

	config.init(
		"cif-diff [options] file1 file2",
		cfg::make_option("help,h", "Display help message"),
		cfg::make_option("version", "Print version"),
		cfg::make_option("verbose,v", "Verbose output"),
		cfg::make_option<std::vector<std::string>>("category", "Limit comparison to this category, default is all categories. Can be specified multiple times"),
		cfg::make_option<int>("max-diff-count", 5, "Maximum number of diff items per category, enter zero (0) for unlimited, default is 5"),
		cfg::make_option("text", "Text based diff (using vimdiff) based on the order of the cif version"),
		cfg::make_option("icase", "Ignore case (vimdiff option)"),
		cfg::make_option("iwhite", "Ignore whitespace (vimdiff option)"),
		cfg::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)")
	);

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().size() != 2)
	{
		std::cerr << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");
	if (config.has("debug"))
		cif::VERBOSE = config.get<int>("debug");

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

	gxrio::ifstream if1{input[0]};
	if (not if1.is_open())
		throw std::runtime_error("Could not open file " + input[0]);

	gxrio::ifstream if2(input[1]);
	if (not if2.is_open())
		throw std::runtime_error("Could not open file " + input[1]);

	cif::file file1{if1};
	cif::file file2{if2};

	if (config.has("text"))
		compareCifsText(file1, file2, fs::path(input[0]), fs::path(input[1]), config.has("icase"), config.has("iwhite"));
	else
		compareCifs(file1.front(), file2.front(), categories, maxDiffCount);

	return 0;
}
