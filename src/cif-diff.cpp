#include "pdb-redo.h"

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

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

// --------------------------------------------------------------------

class templateParser : public cif::SacParser
{
  public:
	templateParser(istream& is)
		: SacParser(is)
	{
	}
	
	virtual void produceDatablock(const string& name)
	{
	}
	
	virtual void produceCategory(const string& name)
	{
	}
	
	virtual void produceRow()
	{
	}
	
	virtual void produceItem(const string& category, const string& item, const string& value)
	{
		string tag = "_" + category + "." + item;
		if (find(mOrder.rbegin(), mOrder.rend(), tag) == mOrder.rend())
			mOrder.push_back(tag);
	}

	vector<string>	mOrder;
};

// --------------------------------------------------------------------

void compareCategories(cif::Category& a, cif::Category& b, size_t maxDiffCount)
{
	using namespace std::placeholders; 
	
//	set<string> tagsA(a.fields()), tagsB(b.fields());
//	
//	if (tagsA != tagsB)
//		cout << "Unequal number of fields" << endl;

	auto& validator = a.getValidator();
	auto catValidator = validator.getValidatorForCategory(a.name());
	if (catValidator == nullptr)
		throw runtime_error("missing cat validator");
	
	typedef function<int(const char*,const char*)> compType;
	vector<tuple<string,compType>> tags;
	auto keys = catValidator->mKeys;
	vector<size_t> keyIx;
	
	for (auto& tag: a.fields())
	{
		auto iv = catValidator->getValidatorForItem(tag);
		if (iv == nullptr)
			throw runtime_error("missing item validator");
		auto tv = iv->mType;
		if (tv == nullptr)
			throw runtime_error("missing type validator");
		tags.push_back(make_tuple(tag, bind(&cif::ValidateType::compare, tv, std::placeholders::_1, std::placeholders::_2)));
		
		auto pred = [tag](const string& s) -> bool { return cif::iequals(tag, s) == 0; };
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
			string tag;
			compType compare;
			
			tie(tag, compare) = tags[kix];

			d = compare(a[tag].c_str(), b[tag].c_str());

			if (d != 0)
				break;
		}
		
		return d < 0;
	};
	
//	vector<cif::Row> rowsA(a.begin(), a.end()), rowsB(b.begin(), b.end());
//	sort(rowsA.begin(), rowsA.end(), rowLess);
//	sort(rowsB.begin(), rowsB.end(), rowLess);
	
	auto ai = a.begin(), bi = b.begin();

	struct Diff
	{
		virtual ~Diff() {}
		
		string key(cif::Row r, vector<string>& keys)
		{
			vector<string> v;
			for (auto k: keys)
				v.push_back(r[k].as<string>());
			return "[" + ba::join(v, ", ") + "]";
		}
		
		virtual void report(vector<string>& keys) = 0;
	};

	struct ExtraADiff : public Diff
	{
		cif::Row A;
		
		ExtraADiff(cif::Row r) : A(r) {}
		
		virtual void report(vector<string>& keys)
		{
			cout << "Extra row in A with key " << key(A, keys) << endl;
		}
	};
	
	struct ExtraBDiff : public Diff
	{
		cif::Row B;
		
		ExtraBDiff(cif::Row r) : B(r) {}
		
		virtual void report(vector<string>& keys)
		{
			cout << "Extra row in B with key " << key(B, keys) << endl;
		}
	};
	
	struct ValueDiff : public Diff
	{
		cif::Row A, B;
		vector<string> missingA, missingB, different;
		
		ValueDiff(cif::Row a, cif::Row b, vector<string>&&  missingA, vector<string>&& missingB, vector<string>&& different)
			: A(a), B(b), missingA(move(missingA)), missingB(move(missingB)), different(move(different)) {}
		
		virtual void report(vector<string>& keys)
		{
			cout << "Differences in rows with key " << key(A, keys) << endl;
			
			for (auto& item: different)
			{
				cout << "    " << item << " (A): '" << A[item] << '\'' << endl
					 << "    " << item << " (B): '" << B[item] << '\'' << endl;
			}
			
			for (auto& item: missingA)
			{
				cout << "    " << item << " (A): <missing>" << endl
					 << "    " << item << " (B): '" << B[item] << '\'' << endl;
			}

			for (auto& item: missingB)
			{
				cout << "    " << item << " (A): '" << A[item] << '\'' << endl
					 << "    " << item << " (B): <missing>" << endl;
			}
		}
	};
	
	vector<Diff*> diffs;

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
		
		vector<string> missingA, missingB, different;
		
		for (auto& tt: tags)
		{
			string tag;
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
		cout << string(cif::get_terminal_width(), '-') << endl
			 << "Differences in values for category " << a.name() << endl
			 << endl;

		for (auto diff: diffs)
		{
			diff->report(keys);
			delete diff;
		}
	
		if (diffs.size() == maxDiffCount)
			cout << "..." << endl;
		cout << endl;
	}
}

void compareCifs(cif::Datablock& dbA, cif::Datablock& dbB, const cif::iset& categories, int maxDiffCount)
{
	vector<string> catA, catB;

	for (auto& cat: dbA)
		catA.push_back(cat.name());
	sort(catA.begin(), catA.end());

	for (auto& cat: dbB)
		catB.push_back(cat.name());
	sort(catB.begin(), catB.end());

	// loop over categories twice, to group output
	// First iteration is to list missing categories.

	vector<string> missingA, missingB;

	auto catA_i = catA.begin(), catB_i = catB.begin();
	
	while (catA_i != catA.end() and catB_i != catB.end())
	{
		string nA = *catA_i;
		ba::to_lower(nA);
		
		string nB = *catB_i;
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
			cout << "Categories missing in A: " << ba::join(missingA, ", ") << endl
				 << endl;
	
		if (not missingB.empty())
			cout << "Categories missing in B: " << ba::join(missingB, ", ") << endl
				 << endl;
	}

	// Second loop, now compare category values
	catA_i = catA.begin(), catB_i = catB.begin();
	
	while (catA_i != catA.end() and catB_i != catB.end())
	{
		string nA = *catA_i;
		ba::to_lower(nA);
		
		string nB = *catB_i;
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
		cerr << "Error creating temp files:  " << strerror(errno) << endl;
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

		vector<string> order;
		a.data().getTagOrder(order);
		b.data().write(out, order);
	}
	
	vector<const char*> nArgv = {
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
			cerr << "Failed to execute vimdiff" << endl;
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
		("category",	po::value<vector<string>>(),	"Limit comparison to this category, default is all categories. Can be specified multiple times")
		("max-diff-count", po::value<int>(),			"Maximum number of diff items per category, enter zero (0) for unlimited, default is 5")
		("text",										"Text based diff (using vimdiff) based on the order of the cif version")
		("icase",										"Ignore case (vimdiff option)")
		("iwhite",										"Ignore whitespace (vimdiff option)");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",	po::value<vector<string>>(),"Input files")
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
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm["input"].as<vector<string>>().size() != 2)
	{
		cerr << visible_options << endl;
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
		for (auto cs: vm["category"].as<vector<string>>())
		{
			for (auto si = ba::make_split_iterator(cs, ba::token_finder(ba::is_any_of(",; "), ba::token_compress_on)); not si.eof(); ++si)
			{
				string cat(si->begin(), si->end());
				ba::to_lower(cat);
				categories.insert(cat);
			}
		}
	}
	
	auto input = vm["input"].as<vector<string>>();
	c::File file1{fs::path(input[0])};
	c::File file2{fs::path(input[1])};

	if (vm.count("text"))
		compareCifsText(file1, file2, vm.count("icase"), vm.count("iwhite"));
	else
		compareCifs(file1.data(), file2.data(), categories, maxDiffCount);

	return 0;	
}

