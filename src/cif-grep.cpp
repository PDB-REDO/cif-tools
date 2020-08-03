#include "pdb-redo.h"

#include <fstream>
#include <functional>
#include <regex>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

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
	grepParser(const string& file, istream& is, const string& pattern, bool quiet, bool printLineNr)
		: SacParser(is), mFile(file), mRx(pattern), mQuiet(quiet), mLineNr(printLineNr)
	{
	}
	
	grepParser(const string& file, istream& is, const string& tag, const string& pattern, bool quiet, bool printLineNr)
		: grepParser(file, is, pattern, quiet, printLineNr)
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
			regex_search(value, mRx))
		{
			++mMatches;
			
			if (not mQuiet)
			{
				if (not mFile.empty())
					cout << mFile << ':';
				if (mLineNr)
					cout << mLineNr << ':';
				cout << value << endl;
			}
		}
	}
	
	string	mFile;
	string	mCat, mItem;
	regex	mRx;
	size_t	mMatches = 0;
	bool	mQuiet, mLineNr;
};

size_t cifGrep(const string& pattern, const string& tag, const string& file, istream& is, bool quiet, bool printLineNr)
{
	size_t result = 0;
	
	if (tag.empty())
	{
		grepParser gp(file, is, pattern, quiet, printLineNr);
		gp.parseFile();
		
		result = gp.getMatches();
	}
	else
	{
		grepParser gp(file, is, tag, pattern, quiet, printLineNr);
		gp.parseFile();
		
		result = gp.getMatches();
	}
	
	return result;
}

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-grep [option...] pattern [file ...]");
	visible_options.add_options()
		("item,i",	po::value<string>(),		"Item tag to scan, default is all item values")
		("help",								"Display help message")
		("version",								"Print version")
		("quiet,q",								"Only print files matching pattern")
		("count,c",								"Only show number of hits")
		("line-number,n",						"Print line numbers")
		("no-filename,h",						"Don't print the filename")
		("with-filename,H",						"Do print the filename")
		("verbose,v",							"Verbose output")
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
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("pattern") == 0)
	{
		cerr << visible_options << endl;
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
	size_t count = 0;

	quiet = quiet or countOnly;

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
		
		if (cif::VERBOSE)
			cerr << "matching only for category: " << cat << " and item " << item << endl;
	}
	
	size_t result = false;
	if (vm.count("input") == 0)
	{
		result = cifGrep(pattern, tag, "stdin", cin, quiet or filenamesOnly, lineNumbers);
		if (doFileNames or (filenamesOnly and result != 0))
			cout << "stdin" << endl;
		if (countOnly)
			cout << result << endl;
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
		
			if (doFileNames)
				noFileNames = false;
			else if (files.size() <= 1)
				noFileNames = true;
		
		for (auto file: filesWithSizes)
		{
			fs::path f;
			size_t size;
			tie(size, f) = file;

			if (not fs::is_regular(f))
				continue;
			
			if (cif::VERBOSE)
				cerr << f << endl;

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
				size_t r = cifGrep(pattern, tag, noFileNames ? ""s : f.filename().string(), in, quiet or filenamesOnly, lineNumbers);

				count += r;

				if (cif::VERBOSE or (countOnly and not noFileNames))
					cout << f << ':' << r << endl;

				if (r > 0)
					result = true;
			}
			catch (const exception& e)
			{
				cerr << endl
						<< "exception for " << f << endl
						<< " => " << e.what() << endl;
			}
		}
	}

	if (noFileNames and countOnly)
		cout << count << endl;

	return result ? 0 : 1;
}

