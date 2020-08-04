#include "pdb-redo.h"

#include <sys/wait.h>

#include <fstream>
#include <chrono>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include <zeep/xml/document.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/Structure.hpp"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;
namespace zx = zeep::xml;

void SwapFields(cif::Row& r, string fld1, string fld2)
{
	string v = r[fld1].as<string>();
	r[fld1] = r[fld2].as<string>();
	r[fld2] = v;
}

// returns <asymId,compId,seqId>
tuple<string,string,string> MapResidue(cif::Datablock& db, string chainID,
	string compID, int seqNum, string iCode)
{
	string asymId, compId, seqId;

	auto r = db["pdbx_poly_seq_scheme"][
		cif::Key("pdb_strand_id") == chainID and
		cif::Key("pdb_mon_id") == compID and
		cif::Key("pdb_seq_num") == seqNum and
		cif::Key("pdb_ins_code") == iCode];
	
	if (not r.empty())
		cif::tie(asymId, compId, seqId) = r.get("asym_id", "mon_id", "seq_id");
	else
	{
		auto r = db["pdbx_nonpoly_scheme"][
			cif::Key("pdb_strand_id") == chainID and
			cif::Key("pdb_mon_id") == compID and
			cif::Key("pdb_seq_num") == seqNum and
			cif::Key("pdb_ins_code") == iCode];

		if (r.empty())
			throw runtime_error("Could not map residue " + chainID + to_string(seqNum) + iCode);
		
		cif::tie(asymId, compId) = r.get("asym_id", "mon_id");
	}

	return make_tuple(asymId, compId, seqId);
}

// returns <chainID,compId,seqNum,iCode>
tuple<char,string,int,char> MapBackResidue(cif::Datablock& db, string asymId,
	string compId, string seqId)
{
	string chainID, iCode;
	int seqNum;

	auto r = db["pdbx_poly_seq_scheme"][
		cif::Key("asym_id") == asymId and
		cif::Key("mon_id") == compId and
		cif::Key("seq_id") == seqId];
	
	if (not r.empty())
		cif::tie(chainID, compId, seqNum, iCode) =
			r.get("pdb_strand_id", "pdb_mon_id", "pdb_seq_num", "pdb_ins_code");
	else
	{
		auto r = db["pdbx_nonpoly_scheme"][
			cif::Key("asym_id") == asymId and
			cif::Key("seq_id") == seqId];

		if (r.empty())
			throw runtime_error("Could not map residue " + asymId + ":" + seqId);
		
		cif::tie(chainID, compId, seqNum, iCode) =
			r.get("pdb_strand_id", "pdb_mon_id", "pdb_seq_num", "pdb_ins_code");
	}

	return make_tuple(chainID[0], compId, seqNum, iCode[0] ? iCode[0] : ' ');
}

// returns <atomId,asymId,monId,seqId>
tuple<string,string,string,string> MapAtom(cif::Datablock& db, string PDB_atomStr)
{
	string atomId = PDB_atomStr.substr(0, 4);		ba::trim(atomId);
	string monId = PDB_atomStr.substr(5, 3);		ba::trim(monId);
	string chainId = PDB_atomStr.substr(9, 1);
	int seqNum = stoi(PDB_atomStr.substr(10, 4));
	string insCode = PDB_atomStr.substr(14, 1);	ba::trim(insCode);
	
	string asymId, seqId;
	
	tie(asymId, monId, seqId) = MapResidue(db, chainId, monId, seqNum, insCode);
	
	return make_tuple(atomId, asymId, monId, seqId);
}


int pr_main(int argc, char* argv[])
{
	int result = 0;
	
	po::options_description visible_options("flipper " + VERSION_STRING + " [options] file]" );
	visible_options.add_options()
		("output,o",	po::value<string>(),	"The output file")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		;

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",		po::value<string>(),	"Input files")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "flipper.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "flipper.conf";
	
	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}
	
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	//
	
	fs::path input = vm["input"].as<string>();
	c::File pdb(input);
	c::Structure structure(pdb);
	
	auto& db = pdb.data();

	const set<string> kDEFY{"TYR", "PHE", "ASP", "GLU"};

	for (auto& poly: structure.polymers())
	{
		for (auto& res: poly)
		{
			string resType = res.compoundID();

			if (not kDEFY.count(resType))
				continue;
			
			size_t chiNr = resType == "GLU" ? 2 : 1;
			if (res.nrOfChis() != chiNr + 1)
				continue;
			
			try
			{
				float angle = res.chi(chiNr);
				if (angle >= -90 and angle <= 90)
					continue;
			}
			catch(const std::exception& e)
			{
				if (cif::VERBOSE)
					std::cerr << e.what() << '\n';
				continue;
			}
			
			//

			cerr << "Need to flip " << res.compoundID() << ' ' << res.asymID() << ':' << res.seqID() << endl;
			
			if (resType == "TYR" or resType == "PHE")
			{
				auto a = res.atomByID("CD1");
				auto b = res.atomByID("CD2");
				structure.swapAtoms(a, b);

				auto c = res.atomByID("CE1");
				auto d = res.atomByID("CE2");
				structure.swapAtoms(c, d);
			}
			else if (resType == "ASP")
			{
				auto a = res.atomByID("OD1");
				auto b = res.atomByID("OD2");
				structure.swapAtoms(a, b);
			}
			else
			{
				auto a = res.atomByID("OE1");
				auto b = res.atomByID("OE2");
				structure.swapAtoms(a, b);
			}
		}
	}

	// done, write out
	db.add_software("flipper", "other", get_version_nr(), get_version_date());

	// structure.sortAtoms();
	// The sort is a 'stable' sort, that means we keep the order unless
	// a change is needed. So ignore sorting on label_atom_id since
	// that will keep the files mostly intact.

	db["atom_site"].sort([](const cif::Row& a, const cif::Row& b) -> int
	{
		int d = 0;

		if (not (a["pdbx_PDB_model_num"].empty() or b["pdbx_PDB_model_num"].empty()))
			d = a["pdbx_PDB_model_num"].as<int>() - b["pdbx_PDB_model_num"].as<int>();

		if (d == 0)
		{
			string asymIDa = a["label_asym_id"].as<string>();
			string asymIDb = b["label_asym_id"].as<string>();

			d = asymIDa.length() - asymIDb.length();
			if (d == 0)
				d = asymIDa.compare(asymIDb);
		}

		if (d == 0)
		{
			int seqIDa = a["label_seq_id"].as<int>();
			int seqIDb = b["label_seq_id"].as<int>();
			d = seqIDa - seqIDb;
		}

		if (d == 0)
		{
			string asymIDa = a["label_asym_id"].as<string>();
			string asymIDb = b["label_asym_id"].as<string>();
			d = asymIDa.compare(asymIDb);
		}


		return d;
	});

	if (vm.count("output"))
	{
		
		
		pdb.save(vm["output"].as<string>());
	}

	return result;
}
