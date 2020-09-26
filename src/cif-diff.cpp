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

#include <sys/wait.h>

#include <fstream>
#include <functional>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Cif2PDB.hpp"
#include "cif++/Structure.hpp"
#include "cif++/CifParser.hpp"
#include "cif++/CifValidator.hpp"
#include "cif++/CifUtils.hpp"

namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

// --------------------------------------------------------------------

class templateParser : public cif::SacParser
{
  public:
	templateParser(std::istream& is)
		: SacParser(is)
	{
	}
	
	virtual void produceDatablock(const std::string& name)
	{
	}
	
	virtual void produceCategory(const std::string& name)
	{
	}
	
	virtual void produceRow()
	{
	}
	
	virtual void produceItem(const std::string& category, const std::string& item, const std::string& value)
	{
		std::string tag = "_" + category + "." + item;
		if (find(mOrder.rbegin(), mOrder.rend(), tag) == mOrder.rend())
			mOrder.push_back(tag);
	}

	std::vector<std::string>	mOrder;
};

// --------------------------------------------------------------------

void compareCategories(cif::Category& a, cif::Category& b, size_t maxDiffCount)
{
	using namespace std::placeholders; 
	
//	set<std::string> tagsA(a.fields()), tagsB(b.fields());
//	
//	if (tagsA != tagsB)
//		std::cout << "Unequal number of fields" << std::endl;

	auto& validator = a.getValidator();
	auto catValidator = validator.getValidatorForCategory(a.name());
	if (catValidator == nullptr)
		throw std::runtime_error("missing cat validator");
	
	typedef std::function<int(const char*,const char*)> compType;
	std::vector<std::tuple<std::string,compType>> tags;
	auto keys = catValidator->mKeys;
	std::vector<size_t> keyIx;
	
	for (auto& tag: a.fields())
	{
		auto iv = catValidator->getValidatorForItem(tag);
		if (iv == nullptr)
			throw std::runtime_error("missing item validator");
		auto tv = iv->mType;
		if (tv == nullptr)
			throw std::runtime_error("missing type validator");
		tags.push_back(std::make_tuple(tag, std::bind(&cif::ValidateType::compare, tv, std::placeholders::_1, std::placeholders::_2)));
		
		auto pred = [tag](const std::string& s) -> bool { return cif::iequals(tag, s) == 0; };
		if (find_if(keys.begin(), keys.end(), pred) == keys.end())
			keyIx.push_back(tags.size() - 1);
	}
	
	a.reorderByIndex();
	b.reorderByIndex();
	
	auto rowLess = [&](const cif::Row& a, const cif::Row& b) -> bool
	{
		int d = 0;

		for (auto kix: keyIx)
		{
			std::string tag;
			compType compare;
			
			tie(tag, compare) = tags[kix];

			d = compare(a[tag].c_str(), b[tag].c_str());

			if (d != 0)
				break;
		}
		
		return d < 0;
	};
	
//	std::vector<cif::Row> rowsA(a.begin(), a.end()), rowsB(b.begin(), b.end());
//	sort(rowsA.begin(), rowsA.end(), rowLess);
//	sort(rowsB.begin(), rowsB.end(), rowLess);
	
	auto ai = a.begin(), bi = b.begin();

	struct Diff
	{
		virtual ~Diff() {}
		
		std::string key(cif::Row r, std::vector<std::string>& keys)
		{
			std::vector<std::string> v;
			for (auto k: keys)
				v.push_back(r[k].as<std::string>());
			return "[" + ba::join(v, ", ") + "]";
		}
		
		virtual void report(std::vector<std::string>& keys) = 0;
	};

	struct ExtraADiff : public Diff
	{
		cif::Row A;
		
		ExtraADiff(cif::Row r) : A(r) {}
		
		virtual void report(std::vector<std::string>& keys)
		{
			std::cout << "Extra row in A with key " << key(A, keys) << std::endl;
		}
	};
	
	struct ExtraBDiff : public Diff
	{
		cif::Row B;
		
		ExtraBDiff(cif::Row r) : B(r) {}
		
		virtual void report(std::vector<std::string>& keys)
		{
			std::cout << "Extra row in B with key " << key(B, keys) << std::endl;
		}
	};
	
	struct ValueDiff : public Diff
	{
		cif::Row A, B;
		std::vector<std::string> missingA, missingB, different;
		
		ValueDiff(cif::Row a, cif::Row b, std::vector<std::string>&&  missingA, std::vector<std::string>&& missingB, std::vector<std::string>&& different)
			: A(a), B(b), missingA(move(missingA)), missingB(move(missingB)), different(move(different)) {}
		
		virtual void report(std::vector<std::string>& keys)
		{
			std::cout << "Differences in rows with key " << key(A, keys) << std::endl;
			
			for (auto& item: different)
			{
				std::cout << "    " << item << " (A): '" << A[item] << '\'' << std::endl
					 << "    " << item << " (B): '" << B[item] << '\'' << std::endl;
			}
			
			for (auto& item: missingA)
			{
				std::cout << "    " << item << " (A): <missing>" << std::endl
					 << "    " << item << " (B): '" << B[item] << '\'' << std::endl;
			}

			for (auto& item: missingB)
			{
				std::cout << "    " << item << " (A): '" << A[item] << '\'' << std::endl
					 << "    " << item << " (B): <missing>" << std::endl;
			}
		}
	};
	
	std::vector<Diff*> diffs;

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
		
		cif::Row ra = *ai, rb = *bi;
		
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
		
		for (auto& tt: tags)
		{
			std::string tag;
			compType compare;
			
			tie(tag, compare) = tt;
			
			// make it an option to compare unapplicable to empty or something
			
			const char* ta = ra[tag].c_str();	if (strcmp(ta, ".") == 0) ta = "";
			const char* tb = rb[tag].c_str();	if (strcmp(tb, ".") == 0) tb = "";
			
			if (compare(ta, tb) != 0)
			{
				if (*ta == 0)
					missingA.push_back(tag);
				else if (*tb == 0)
					missingB.push_back(tag);
				else
					different.push_back(tag);
			}
		}
		
		++ai;
		++bi;
		
		if (not missingA.empty() or not missingB.empty() or not different.empty())
			diffs.push_back(new ValueDiff{ ra, rb, move(missingA), move(missingB), move(different) });
	}

	if (not diffs.empty())
	{
		std::cout << std::string(cif::get_terminal_width(), '-') << std::endl
			 << "Differences in values for category " << a.name() << std::endl
			 << std::endl;

		for (auto diff: diffs)
		{
			diff->report(keys);
			delete diff;
		}
	
		if (diffs.size() == maxDiffCount)
			std::cout << "..." << std::endl;
		std::cout << std::endl;
	}
}

void compareCifs(cif::Datablock& dbA, cif::Datablock& dbB, const cif::iset& categories, int maxDiffCount)
{
	std::vector<std::string> catA, catB;

	for (auto& cat: dbA)
		catA.push_back(cat.name());
	sort(catA.begin(), catA.end());

	for (auto& cat: dbB)
		catB.push_back(cat.name());
	sort(catB.begin(), catB.end());

	// loop over categories twice, to group output
	// First iteration is to list missing categories.

	std::vector<std::string> missingA, missingB;

	auto catA_i = catA.begin(), catB_i = catB.begin();
	
	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		ba::to_lower(nA);
		
		std::string nB = *catB_i;
		ba::to_lower(nB);
		
		int d = nA.compare(nB);
		if (d > 0)
		{
			auto& cat = dbB[*catB_i++];
			
			if (not cat.empty())
				missingA.push_back(cat.name());
		}
		else if (d < 0)
		{
			auto& cat = dbA[*catA_i++];
			
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
			std::cout << "Categories missing in A: " << ba::join(missingA, ", ") << std::endl
				 << std::endl;
	
		if (not missingB.empty())
			std::cout << "Categories missing in B: " << ba::join(missingB, ", ") << std::endl
				 << std::endl;
	}

	// Second loop, now compare category values
	catA_i = catA.begin(), catB_i = catB.begin();
	
	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		ba::to_lower(nA);
		
		std::string nB = *catB_i;
		ba::to_lower(nB);
		
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

void compareCifsText(c::File& a, c::File& b, bool icase, bool iwhite)
{
	// temp files for vimdiff

	char generated[] = "/tmp/pdb2cif-diff-B-XXXXXX.cif", original[] = "/tmp/pdb2cif-diff-A-XXXXXX.cif";
	int fd[2];
	
	if ((fd[0] = mkstemps(generated, 4)) < 0 or (fd[1] = mkstemps(original, 4)) < 0)
	{
		std::cerr << "Error creating temp files:  " << strerror(errno) << std::endl;
		exit(1);
	}

	io::file_descriptor_sink gen(fd[0], io::close_handle), orig(fd[1], io::close_handle);

	{
		io::filtering_stream<io::output> out;
		out.push(orig);

		a.data().write(out);
	}

	// Next the converted cif file
	
	{
		io::filtering_stream<io::output> out;
		out.push(gen);

		std::vector<std::string> order;
		a.data().getTagOrder(order);
		b.data().write(out, order);
	}
	
	std::vector<const char*> nArgv = {
		"/usr/bin/vimdiff"
	};
	
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

	nArgv.push_back(original);
	nArgv.push_back(generated);
	nArgv.push_back(nullptr);
	
	int pid = fork();
	
	if (pid <= 0)
	{
		if (execv(nArgv[0], const_cast<char*const*>(nArgv.data())) < 0)
			std::cerr << "Failed to execute vimdiff" << std::endl;
		exit(1);
	}
	
	int status;
	waitpid(pid, &status, 0);
	
	if (WIFEXITED(status))
	{
		unlink(generated);
		unlink(original);
	}	
}

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-diff " + VERSION_STRING + " options file1 file2");
	visible_options.add_options()
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		("category",	po::value<std::vector<std::string>>(),	"Limit comparison to this category, default is all categories. Can be specified multiple times")
		("max-diff-count", po::value<int>(),			"Maximum number of diff items per category, enter zero (0) for unlimited, default is 5")
		("text",										"Text based diff (using vimdiff) based on the order of the cif version")
		("icase",										"Ignore case (vimdiff option)")
		("iwhite",										"Ignore whitespace (vimdiff option)");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",	po::value<std::vector<std::string>>(),"Input files")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm["input"].as<std::vector<std::string>>().size() != 2)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	int maxDiffCount = 5;
	if (vm.count("max-diff-count"))
		maxDiffCount = vm["max-diff-count"].as<int>();
	
	cif::iset categories;
	if (vm.count("category"))
	{
		for (auto cs: vm["category"].as<std::vector<std::string>>())
		{
			for (auto si = ba::make_split_iterator(cs, ba::token_finder(ba::is_any_of(",; "), ba::token_compress_on)); not si.eof(); ++si)
			{
				std::string cat(si->begin(), si->end());
				ba::to_lower(cat);
				categories.insert(cat);
			}
		}
	}
	
	auto input = vm["input"].as<std::vector<std::string>>();
	c::File file1{fs::path(input[0])};
	c::File file2{fs::path(input[1])};

	if (vm.count("text"))
		compareCifsText(file1, file2, vm.count("icase"), vm.count("iwhite"));
	else
		compareCifs(file1.data(), file2.data(), categories, maxDiffCount);

	return 0;	
}

