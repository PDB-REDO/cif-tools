#include "pdb-redo.h"

#include <fstream>
#include <functional>
#include <regex>
#include <atomic>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/thread.hpp>

#include "cif++/Cif++.h"
#include "cif++/Structure.h"
#include "cif++/CifParser.h"
#include "cif++/CifUtils.h"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

class grepParser : public cif::SacParser
{
  public:
	grepParser(const string& file, istream& is, const string& pattern, bool quiet)
		: SacParser(is), mFile(file), mRx(pattern), mQuiet(quiet)
	{
	}
	
	grepParser(const string& file, istream& is, const string& tag, const string& pattern, bool quiet)
		: SacParser(is), mFile(file), mRx(pattern), mQuiet(quiet)
	{
		tie(mCat, mItem) = cif::splitTagName(tag);
	}
	
	size_t getMatches() const			{ return mMatches; }
	
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
		if ((mCat.empty() or cif::iequals(category, mCat)) and
			(mItem.empty() or cif::iequals(item, mItem)) and
			regex_match(value, mRx))
		{
			++mMatches;
			
			if (not mQuiet)
			{
				if (not mFile.empty())
					cout << mFile << ':';
				cout << mLineNr << '\t' << value << endl;
			}
		}
	}
	
	string	mFile;
	string	mCat, mItem;
	regex	mRx;
	size_t	mMatches = 0;
	bool	mQuiet;
};

bool cifGrep(const string& pattern, const string& tag, const string& file, istream& is, bool quiet)
{
	bool result = false;
	
	if (tag.empty())
	{
		grepParser gp(file, is, pattern, quiet);
		gp.parseFile();
		
		result = gp.getMatches() > 0;
	}
	else
	{
		grepParser gp(file, is, tag, pattern, quiet);
		gp.parseFile();
		
		result = gp.getMatches() > 0;
	}
	
	return result;
}

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-grep [option...] pattern [file ...]");
	visible_options.add_options()
		("item,i",	po::value<string>(),		"Item tag to scan, default is all item values")
		("help,h",								"Display help message")
		("version",								"Print version")
		("quiet",								"Only print files matching pattern")
		("verbose,v",							"Verbose output")
		("progress",							"Show progress bar")
		("files-with-matches,l",				"Print only names of files containing matches")
		("recursive,r",							"Search recursively");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("pattern",	po::value<string>(),		"Pattern")
		("input",	po::value<vector<string>>(),"Input files")
		("debug,d",	po::value<int>(),			"Debug level (for even more verbose output)");

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
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("pattern") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();

	bool quiet = vm.count("quiet") > 0;
	bool filenamesOnly = vm.count("files-with-matches") > 0;
	bool progress = VERBOSE == 0 and vm.count("progress") > 0;

	string pattern = vm["pattern"].as<string>();
	string tag;
	if (vm.count("item"))
	{
		tag = vm["item"].as<string>();

		string cat, item;
		tie(cat, item) = cif::splitTagName(tag);
		
		if (cat.empty())
			throw runtime_error("Invalid category in tag: '" + cat + '\'');
		
		if (item.empty())
			throw runtime_error("Invalid item: '" + item + '\''); 
		
		if (VERBOSE)
			cerr << "matching only for category: " << cat << " and item " << item << endl;
	}
	
	bool result = false;
	if (vm.count("input") == 0)
	{
		result = cifGrep(pattern, tag, "stdin", cin, quiet or filenamesOnly);
		if (filenamesOnly and result)
			cout << "stdin" << endl;
	}
	else
	{
		auto files = vm["input"].as<vector<string>>();
		
		if (vm.count("recursive"))
		{
			vector<string> expanded;
			for (auto file: files)
			{
				if (fs::is_directory(file))
				{
					for (auto i = fs::recursive_directory_iterator(file);
						i != fs::recursive_directory_iterator(); ++i)
					{
						fs::path p = i->path();
						if (fs::is_regular(p))
							expanded.push_back(p.string());
					}
				}
				else
					expanded.push_back(file);
			}
			
			files = expanded;
		}
		
		vector<tuple<size_t,string>> filesWithSizes;
		size_t totalSize = 0;
		
		transform(files.begin(), files.end(), back_inserter(filesWithSizes),
			[&totalSize](const string& f) -> tuple<size_t, string>
			{
				size_t size = fs::file_size(f);
				totalSize += size;
				return make_tuple(size, f);
			});
		
//		sort(filesWithSizes.begin(), filesWithSizes.end(),
//			[](auto& a, auto& b) -> bool
//			{
//				return get<0>(a) > get<0>(b);
//			});
		
		unique_ptr<cif::Progress> p;
		if (progress)
			p.reset(new cif::Progress(totalSize, "searching"));
		
		vector<fs::path> matches;
		
		boost::thread_group t;
		size_t N = boost::thread::hardware_concurrency();
		atomic<size_t> next(0);
		
		for (size_t i = 0; i < N; ++i)
			t.create_thread([&]()
			{
				for (;;)
				{
					size_t i = next++;
					
					if (i >= filesWithSizes.size())
						break;

					size_t size;
					fs::path f;
					tie(size, f) = filesWithSizes[i];
					
					if (not fs::is_regular(f))
						continue;
					
					if (VERBOSE)
						cerr << f << endl;

					if (p)
						p->message(f.leaf().string());
					
					ifstream infile(f.c_str(), ios_base::in | ios_base::binary);
					if (not infile.is_open())
						throw runtime_error("Could not open file " + f.string());
			
					io::filtering_stream<io::input> in;
				
					if (f.extension() == ".bz2")
						in.push(io::bzip2_decompressor());
					else if (f.extension() == ".gz")
						in.push(io::gzip_decompressor());
					
					in.push(infile);
			
					try
					{
						if (cifGrep(pattern, tag, f.filename().string(), in, quiet or filenamesOnly))
						{
							matches.push_back(f);
							if (filenamesOnly)
								cout << f << endl;
							else if (VERBOSE)
								cout << "Matching file: " << f << endl;
							result = true;
						}
					}
					catch (const exception& e)
					{
						cerr << endl
							 << "exception for " << f << endl
							 << " => " << e.what() << endl;
					}
						
					if (p)
						p->consumed(size);
				}
			});
		
		t.join_all();
	}

	return result ? 0 : 1;
}

