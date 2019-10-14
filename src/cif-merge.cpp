#include "pdb-redo.h"

#include <sys/wait.h>

#include <fstream>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include "cif++/Cif++.h"
#include "cif++/Cif2PDB.h"
#include "cif++/Structure.h"
#include "cif++/CifParser.h"
#include "cif++/CifValidator.h"
#include "cif++/CifUtils.h"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

using cif::iequals;

// --------------------------------------------------------------------

void updateEntryID(cif::File& target, const string& entryID)
{
	auto& db = target.firstDatablock();
	
	if (db.getName() != entryID)
		db.setName(entryID);

	for (auto r: db["entry"])
		r["id"] = entryID;
}

void transplant(cif::File& target, cif::File& donor)
{
	auto& dbt = target.firstDatablock();
	auto& dbd = donor.firstDatablock();

	for (auto c: { "struct", "struct_keywords", "audit_author", "citation", "citation_author", "diffrn", "diffrn_radiation", "diffrn_radiation_wavelength" })
	{
		auto cd = dbd.get(c);
		if (cd == nullptr or cd->empty())
			continue;

		auto ct = dbt.get(c);
		if (ct == nullptr)
		{
			dbt.emplace(c);
			ct = dbt.get(c);
		}
		
		ct->clear();
		
		for (auto r: *cd)
			ct->emplace(r);
	}
	
	string exptlMethod;
	auto dExplt = dbd["exptl"].find(cif::Key("entry_id") == dbd.getName());
	if (dExplt.size() != 1)
		throw runtime_error("Invalid number of exptl records in donor file, should be exactly one this version of cif-merge");
	cif::tie(exptlMethod) = dExplt.front().get("method");
	auto tExplt = dbt["exptl"].find(cif::Key("entry_id") == dbt.getName());
	if (tExplt.empty())
	{
		auto c = dbt.emplace("exptl");
		get<0>(c)->emplace({
			{ "entry_id", dbt.getName() },
			{ "method", exptlMethod }
		});
	}
	else
		tExplt.front()["method"] = exptlMethod;
	
	// create a mapping for the entity_ids in both files

	const map<string,const char*> kSrcMap{
		{ "man", "entity_src_gen" },
		{ "nat", "entity_src_nat" },
		{ "syn", "pdbx_entity_src_syn" }
	};
	
	map<string,string> d2tEntityIds;
	auto& targetEntity = dbt["entity"];
	
	for (auto r : targetEntity)
	{
		string id, type, dEntityID;
		cif::tie(id, type) = r.get("id", "type");
		
		if (iequals(type, "polymer"))
		{
			auto t = dbt["entity_poly"][cif::Key("entity_id") == id];
			
			string polyType, seq;
			cif::tie(polyType, seq) = t.get("type", "pdbx_seq_one_letter_code");
			
			auto d = dbd["entity_poly"][cif::Key("type") == polyType and cif::Key("pdbx_seq_one_letter_code") == seq];
			
			if (d.empty())
			{
				if (cif::VERBOSE)
					cerr << "Cannot map entity " << id << " in target file to an entity in the donor" << endl;
				continue;
			}
			
			dEntityID = d["entity_id"].as<string>();
			
			// copy over refseq
			
			auto sr = dbd["struct_ref"][cif::Key("entity_id") == dEntityID];
			if (not sr.empty())
			{
				sr["entity_id"] = id;
				dbt["struct_ref"].emplace(sr);
				
				string refID = sr["id"].as<string>();
				
				for (auto r: dbd["struct_ref_seq"].find(cif::Key("ref_id") == refID))
					dbt["struct_ref_seq"].emplace(r);
		}
		}
		else if (iequals(type, "non-polymer"))
		{
			auto t = dbt["pdbx_entity_nonpoly"][cif::Key("entity_id") == id];
			
			string compID;
			cif::tie(compID) = t.get("comp_id");
			
			auto d = dbd["pdbx_entity_nonpoly"][cif::Key("comp_id") == compID];
			
			if (d.empty())
			{
				if (cif::VERBOSE)
					cerr << "Cannot map entity " << id << " in target file to an entity in the donor" << endl;
				continue;
			}
			
			cif::tie(dEntityID) = d.get("entity_id");
		}
		else if (iequals(type, "water"))
		{
			cif::tie(dEntityID) = dbd["entity"][cif::Key("type") == type].get("id");
		}
		else if (cif::VERBOSE)
			cerr << "Unsupported entity type: " << type << endl;
		
		if (dEntityID.empty())
			continue;

		string srcMethod, description, weight;
		cif::tie(srcMethod, description, weight) =
			dbd["entity"][cif::Key("id") == dEntityID].get("src_method", "pdbx_description", "formula_weight");
		
		r["src_method"] = srcMethod;
		r["pdbx_description"] = description;
		r["formula_weight"] = weight;
		
		if (kSrcMap.count(srcMethod))
		{
			string srcRec = kSrcMap.at(srcMethod);

			auto d = dbd[srcRec][cif::Key("entity_id") == dEntityID];
			if (not d.empty())
			{
				d["entity_id"] = id;
				dbt.emplace(srcRec);
				dbt[srcRec].emplace(d);
			}
		}
	}
}

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-merge " + VERSION + " [options] inputFile donorFile ");
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("input,i",		po::value<string>(),	"Modified PDB file")
		("output,o",	po::value<string>(),	"Output file, default is stdout (terminal)")
		("donor",		po::value<string>(),	"CIF file (or PDB ID for this file) containing the data to collect data from")
		("dict",		po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
		;
	

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("donor", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm.count("donor") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// Load dict, if any
	
	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());

	// Read input file
	mmcif::File cf{vm["input"].as<string>()};
	
	// Read donor file
	mmcif::File df{vm["donor"].as<string>()};

	updateEntryID(cf.file(), df.data().getName());
	transplant(cf.file(), df.file());
	
	if (vm.count("output"))
		cf.save(vm["output"].as<string>());
	else
		cf.file().save(cout);

	return 0;	
}
