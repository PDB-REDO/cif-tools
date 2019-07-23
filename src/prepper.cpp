#include "pdb-redo.h"

#include <sys/wait.h>

#include <fstream>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include <zeep/xml/document.hpp>

#include "cif++/Cif++.h"
#include "cif++/Compound.h"
#include "cif++/Structure.h"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;
namespace zx = zeep::xml;

const set<string> kStripperResidues = {
	"GLY", "ALA", "VAL", "THR", "SER", "CYS", "TRP", "PHE", "TYR", "HIS", "ILE",
	"LEU", "ASP", "ASN", "GLU", "GLN", "MET", "ARG", "LYS", "PRO", "MSE"
};

const set<string> kBackBone = {
	"N", "CA", "C", "O", "OXT"
};

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

int DropLink(cif::Datablock& db, 
	string fromAtomId, string fromAsymId, string fromCompId, string fromSeqId,
	string toAtomId, string toAsymId, string toCompId, string toSeqId)
{
	int result = 0;
	
	for (auto r: db["struct_conn"].find(
		(
			cif::Key("ptnr1_label_asym_id") == fromAsymId and
			cif::Key("ptnr1_label_comp_id") == fromCompId and
			cif::Key("ptnr1_label_seq_id") == fromSeqId and
			cif::Key("ptnr1_label_atom_id") == fromAtomId and
			
			cif::Key("ptnr2_label_asym_id") == toAsymId and
			cif::Key("ptnr2_label_comp_id") == toCompId and
			cif::Key("ptnr2_label_seq_id") == toSeqId and
			cif::Key("ptnr2_label_atom_id") == toAtomId
		) or (
			cif::Key("ptnr1_label_asym_id") == toAsymId and
			cif::Key("ptnr1_label_comp_id") == toCompId and
			cif::Key("ptnr1_label_seq_id") == toSeqId and
			cif::Key("ptnr1_label_atom_id") == toAtomId and
			
			cif::Key("ptnr2_label_asym_id") == fromAsymId and
			cif::Key("ptnr2_label_comp_id") == fromCompId and
			cif::Key("ptnr2_label_seq_id") == fromSeqId and
			cif::Key("ptnr2_label_atom_id") == fromAtomId
		)))
	{
		if (VERBOSE)
			cerr << "Stripping link " << r["id"].as<string>() << endl;
		
		db["struct_conn"].erase(r);
		++result;
	}
	
	return result;
}

tuple<uint32_t,uint32_t> HandlePDBCare(cif::Datablock& db, fs::path pdbCareFile, cif::Datablock& dat)
{
	fs::ifstream file(pdbCareFile);
	if (not file.is_open())
		throw runtime_error("Could not open pdb-care file " + pdbCareFile.string());
	
	zx::document care(file);
	
	int deletedLinkCount = 0, addedLinkCount = 0;
	
	// non_reducingend_links_to_protein handling
	for (auto e: care.find("//non_reducingend_links_to_protein/link"))
	{
		string fromAtomId, fromAsymId, fromCompId, fromSeqId;
		tie(fromAtomId, fromAsymId, fromCompId, fromSeqId)
			= MapAtom(db, e->get_attribute("from"));

		string toAtomId, toAsymId, toCompId, toSeqId;
		tie(toAtomId, toAsymId, toCompId, toSeqId)
			= MapAtom(db, e->get_attribute("to"));
		
		deletedLinkCount += DropLink(db, fromAtomId, fromAsymId, fromCompId, fromSeqId,
			toAtomId, toAsymId, toCompId, toSeqId);
	}
	
	// delete_link
	for (auto e: care.find("//delete_link"))
	{
		string fromAtomId, fromAsymId, fromCompId, fromSeqId;
		tie(fromAtomId, fromAsymId, fromCompId, fromSeqId)
			= MapAtom(db, e->get_attribute("from"));

		string toAtomId, toAsymId, toCompId, toSeqId;
		tie(toAtomId, toAsymId, toCompId, toSeqId)
			= MapAtom(db, e->get_attribute("to"));
		
		deletedLinkCount += DropLink(db, fromAtomId, fromAsymId, fromCompId, fromSeqId,
			toAtomId, toAsymId, toCompId, toSeqId);
	}

	// create_link
	int linkID = 1;
	for (auto e: care.find("//create_link"))
	{
		string fromAtomId, fromAsymId, fromCompId, fromSeqId;
		tie(fromAtomId, fromAsymId, fromCompId, fromSeqId)
			= MapAtom(db, e->get_attribute("from"));

		string toAtomId, toAsymId, toCompId, toSeqId;
		tie(toAtomId, toAsymId, toCompId, toSeqId)
			= MapAtom(db, e->get_attribute("to"));
		
		db["struct_conn"].emplace({
			{ "id", "covale_s" + to_string(linkID++) },
			{ "conn_type_id", "covale" },
			{ "ptnr1_label_asym_id", fromAsymId },
			{ "ptnr1_label_comp_id", fromCompId },
			{ "ptnr1_label_seq_id", fromSeqId.empty() ? "." : fromSeqId },
			{ "ptnr1_label_atom_id", fromAtomId },
			
			{ "ptnr2_label_asym_id", toAsymId },
			{ "ptnr2_label_comp_id", toCompId },
			{ "ptnr2_label_seq_id", toSeqId.empty() ? "." : toSeqId },
			{ "ptnr2_label_atom_id", toAtomId }
		});
		
		++addedLinkCount;
	}
	
	set<string> oldNames, newNames;
	
	for (auto e: care.find("//issue_remedies/node[@issues > 0]"))
	{
		auto rename = e->find_first("remedy/alternative/rename_residue");
		if (not rename)
			continue;
		
		string pdbResname = e->get_attribute("pdb_resname");
		
		string resname = pdbResname.substr(0, 3);
		string chainID = pdbResname.substr(4, 1);
		int seqNr = stoi(pdbResname.substr(5, 4));
		string iCode = pdbResname.substr(9, 1);
		
		if (iCode == " ")
			iCode.clear();
		
		// see if we need to skip this pair
		string newName = rename->get_attribute("new_name");
		if (not dat["carb_rename_ignore"].find(cif::Key("from_comp_id") == resname and cif::Key("to_comp_id") == newName).empty())
		{
			if (VERBOSE)
				cerr << "Skipping pdb-care rename suggestion from " << resname << " to " << newName << endl;
			continue;
		}

		// Now map back this PDB residue to a mmCIF residue
		
		auto ppss = db["pdbx_nonpoly_scheme"].find(
			cif::Key("pdb_strand_id") == chainID and
			cif::Key("pdb_mon_id") == resname and
			cif::Key("pdb_seq_num") == seqNr and
			cif::Key("pdb_ins_code") == iCode);
		
		if (ppss.empty())
		{
			cerr << "Could not rename '" << pdbResname << "' since it was not found" << endl;
			continue;
		}
		
		assert(ppss.size() == 1);
		
		string asymId, entityId, monId;
		cif::tie(asymId, entityId, monId) =
			ppss.front().get("asym_id", "entity_id", "mon_id");
			
		// OK, now we can rename the residue
		
		ppss.front()["mon_id"] = newName;
		
		// And additional records
		for (auto r: db["atom_site"].find(
			cif::Key("label_comp_id") == monId and
			cif::Key("label_asym_id") == asymId and
			cif::Key("label_entity_id") == entityId))
		{
			r["label_comp_id"] = newName;
			
			// is this wise?
			r["auth_comp_id"] = newName;
		}
		
		if (VERBOSE)
			cout << "Renamed " << pdbResname << " to " << newName << endl;
		
		oldNames.insert(resname);
		newNames.insert(newName);
	}
	
	// remove chem_comp for oldNames, if not referenced by any atom_site
	for (auto cmp: oldNames)
	{
		if (db["atom_site"].find(cif::Key("label_comp_id") == cmp).empty())
		{
			if (VERBOSE)
				cout << "Removing _chem_comp " << cmp << endl;

			db["chem_comp"].erase(cif::Key("id") == cmp);
		}
	}

	for (auto cmp: newNames)
	{
		if (not db["chem_comp"].find(cif::Key("id") == cmp).empty())
			continue;

		if (VERBOSE)
			cout << "Adding _chem_comp " << cmp << endl;
		
		auto compound = c::Compound::create(cmp);
		if (compound == nullptr)
		{
			if (VERBOSE)
				cerr << "Unknown compound " << cmp << endl;

			db["chem_comp"].emplace({
				{ "id", cmp },
				{ "mon_nstd_flag", "." },
				{ "type", "saccharide" }
			});
			
			continue;
		}
		
		string formula = compound->formula();
		string name = compound->name();
		string type = compound->type();
		
		if (type.empty())
			type = "saccharide";

		db["chem_comp"].emplace({
			{ "id", cmp },
			{ "name", name },
			{ "formula", formula },
			{ "mon_nstd_flag", "." },
			{ "type", type }
		});
	}
	
	db["chem_comp"].reorderByIndex();
	
//	cout << deletedLinkCount << " carbohydrate atoms in total were deleted." << endl
//		 << 
//      WRITE(6,*) LNKKIL, ' carbohydrate LINKs in total will be deleted.'
//      WRITE(6,*) NALINK, ' carbohydrate LINKs in total will be added.'
//      WRITE(6,*) RENCNT, ' carbohydrate residues will be renamed.'
//C      PRINT*, LINLNK

	return make_tuple(deletedLinkCount, addedLinkCount);
}

void FixMatrix(cif::Datablock& db)
{
	const array<float,12> kIdentity({
		1, 0, 0,	// matrix1
		0, 1, 0,	// matrix2
		0, 0, 1,	// matrix3
		0, 0, 0		// vector
	});
	
	bool hasIdentiy = false;
	array<float,12> mat;

	auto& structNcsOper = db["struct_ncs_oper"];
	
	for (auto r: structNcsOper)
	{
		cif::tie(mat[0], mat[1], mat[2], mat[3], mat[4], mat[5],
				 mat[6], mat[7], mat[8], mat[9], mat[10], mat[11])
			= r.get("matrix[1][1]", "matrix[1][2]", "matrix[1][3]",
					"matrix[2][1]", "matrix[2][2]", "matrix[2][3]",
					"matrix[3][1]", "matrix[3][2]", "matrix[3][3]",
					"vector[1]", "vector[2]", "vector[3]");
		
		// force given flag on identity matrix
		if (mat == kIdentity)
		{
			r["code"] = "given";
			hasIdentiy = true;
			break;
		}
	}
	
	if (not hasIdentiy and not structNcsOper.empty())
	{
		// Need to add the identity matrix
		
		structNcsOper.emplace({
			{ "id", 0 },
			{ "code", "given" },
			{ "matrix[1][1]", kIdentity[0] }, 
			{ "matrix[1][2]", kIdentity[1] }, 
			{ "matrix[1][3]", kIdentity[2] }, 
			{ "matrix[2][1]", kIdentity[3] }, 
			{ "matrix[2][2]", kIdentity[4] }, 
			{ "matrix[2][3]", kIdentity[5] }, 
			{ "matrix[3][1]", kIdentity[6] }, 
			{ "matrix[3][2]", kIdentity[7] }, 
			{ "matrix[3][3]", kIdentity[8] }, 
			{ "vector[1]", 0 }, 
			{ "vector[2]", 0 }, 
			{ "vector[3]", 0 } 
		});
	}
	
	// strip identity scale, if present (is this really necessary?)
	for (auto r: db["atom_sites"])
	{
		cif::tie(mat[0], mat[1], mat[2], mat[3], mat[4], mat[5],
				 mat[6], mat[7], mat[8], mat[9], mat[10], mat[11])
			= r.get("fract_transf_matrix[1][1]", "fract_transf_matrix[1][2]", "fract_transf_matrix[1][3]",
					"fract_transf_matrix[2][1]", "fract_transf_matrix[2][2]", "fract_transf_matrix[2][3]",
					"fract_transf_matrix[3][1]", "fract_transf_matrix[3][2]", "fract_transf_matrix[3][3]",
					"fract_transf_vector[1]", "fract_transf_vector[2]", "fract_transf_vector[3]");
		
		if (mat == kIdentity)
		{
			db["atom_sites"].erase(r);
			break;
		}
	}
}

int pr_main(int argc, char* argv[])
{
	int result = 0;
	
	po::options_description visible_options("stripper " + VERSION + " options file]" );
	visible_options.add_options()
		("output,o",	po::value<string>(),	"The output file, default is stdout")
		("pdb-care",	po::value<string>(),	"pdb-care file")
		("help,h",								"Display help message")
		("version",								"Print version")
		("validate",							"Validate the content of the file using the default or specified dictionary")
		("verbose,v",							"Verbose output")
		("server",								"Server mode")
		("pdb-redo-data", po::value<string>(),	"The PDB-REDO dat file" /*, default is the built in one"*/)
		("dict",		po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
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

	fs::path configFile = "prepper.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "prepper.conf";
	
	if (fs::exists(configFile))
	{
		fs::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}
	
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	fs::path pdbRedoDataFile = "pdb-redo-data.cif";
	if (vm.count("pdb-redo-data"))
		pdbRedoDataFile = vm["pdb-redo-data"].as<string>();
	fs::ifstream pdbRedoData(pdbRedoDataFile);
	if (not pdbRedoData.is_open())
		throw runtime_error("Could not open pdb-redo-data file");
	cif::File pdbRedoDataCif(pdbRedoData);
	cif::Datablock& dat = pdbRedoDataCif.firstDatablock();

	// Load dict, if any
	
	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();
	
	bool serverMode = vm.count("server") > 0;
	
	fs::path input = vm["input"].as<string>();
	c::File pdb(input);
	
	auto& db = pdb.data();

	auto entityForAsym = [&db](const string& asymId)
	{
		string entityId;
		for (auto& r: db["struct_asym"].find(cif::Key("id") == asymId))
		{
			entityId = r["entity_id"].as<string>();
			break;
		}
		return entityId;
	};

	
	// Ready to do some work

	uint32_t numOfOccupanciesReset = 0,
		numOfAtomsDeleted = 0,
		numOfLinksDeleted = 0,
		numOfLinksAdded = 0,
		numOfLinksFixed = 0,
		numOfAnisosDeleted = 0,
		numOfBFactorsReset = 0,
		numOfRenamedWaters = 0;

	FixMatrix(db);
	
	if (vm.count("pdb-care"))
		tie(numOfLinksDeleted, numOfLinksAdded) = HandlePDBCare(db, vm["pdb-care"].as<string>(), dat);
	
	// handle atoms with zero occupancy
	for (auto a: db["atom_site"].find(cif::Key("occupancy") == 0.0))
	{
		string id, compId, atomId, altId;
		cif::tie(id, compId, atomId, altId) = a.get("id", "label_comp_id", "label_atom_id", "label_alt_id");
		
		if (compId == "HOH")
		{
			if (VERBOSE)
				cerr << "Deleted zero occupancy water atom: " << a["id"] << endl;

			db["atom_site"].erase(a);
			++numOfAtomsDeleted;
		}
		else
		{
			auto compound = c::Compound::create(compId);
			if (compound and (cif::iequals(compound->group(), "peptide") or cif::iequals(compound->group(), "P-peptide")))
			{
				if (VERBOSE)
					cerr << "Changed occupancy to 1.00 for (hetero) atom: " << id << endl;
				a["occupancy"] = 1.0;
	
				++numOfOccupanciesReset;
			}
			else
			{
				if (VERBOSE)
					cerr << "Changed occupancy to 0.01 for (hetero) atom: " << id << endl;
				a["occupancy"] = 0.01;
	
				++numOfOccupanciesReset;
			}
		}
	}
	
	// handle atoms with negative occupancy
	for (auto a: db["atom_site"].find(cif::Key("occupancy") < 0.0))
	{
		string id, compId, atomId, altId;
		cif::tie(id, compId, atomId, altId) = a.get("id", "label_comp_id", "label_atom_id", "label_alt_id");
		
		if (VERBOSE)
			cerr << "Deleted zero occupancy (hetero) atom: " << a["id"] << endl;

		db["atom_site"].erase(a);
		++numOfAtomsDeleted;
	}
	
	// remove dubious links
	for (auto l: dat["link_remove"])
	{
		string atom[2], residue[2];
		cif::tie(atom[0], residue[0], atom[1], residue[1]) = l.get("atom_1", "residue_1", "atom_2", "residue_2");
		
		auto ll = db["struct_conn"].find(
			(
				cif::Key("ptnr1_label_comp_id") == residue[0] and
				cif::Key("ptnr1_label_atom_id") == atom[0] and
				cif::Key("ptnr2_label_comp_id") == residue[1] and
				cif::Key("ptnr2_label_atom_id") == atom[1]
			) or (
				cif::Key("ptnr2_label_comp_id") == residue[0] and
				cif::Key("ptnr2_label_atom_id") == atom[0] and
				cif::Key("ptnr1_label_comp_id") == residue[1] and
				cif::Key("ptnr1_label_atom_id") == atom[1]
			)
		);
		
		for (auto lr: ll)
		{
			if (VERBOSE)
				cerr << "Deleted dubious LINK " << lr["id"] << endl;

			db["struct_conn"].erase(lr);

			++numOfLinksDeleted;
		}
	}
	
	// remove rediculously long LINKs
	for (auto lr: db["struct_conn"].find(cif::Key("pdbx_dist_value") > 5.0))
	{
		if (VERBOSE)
			cerr << "Deleted extremely long LINK: " << lr["id"] << endl;

		db["struct_conn"].erase(lr);

		++numOfLinksDeleted;
	}
	
	// Flip ASN if glycosylated at OD1 (LINK part)
	for (auto lr: db["struct_conn"].find(
			(
				cif::Key("ptnr1_label_comp_id") == "ASN" and
				cif::Key("ptnr1_label_atom_id") == "OD1" and
				cif::Key("ptnr2_label_comp_id") == "NAG" and
				cif::Key("ptnr2_label_atom_id") == "C1"
			) or (
				cif::Key("ptnr2_label_comp_id") == "ASN" and
				cif::Key("ptnr2_label_atom_id") == "OD1" and
				cif::Key("ptnr1_label_comp_id") == "NAG" and
				cif::Key("ptnr1_label_atom_id") == "C1"
			)
		))
	{
		string asymId, seqId;
		
		if (VERBOSE)
			cerr << "Flipping at OD1 glycosylated ASN " << lr["id"] << endl;
		
		if (lr["ptnr1_label_atom_id"].as<string>() == "OD1")
		{
			cif::tie(asymId, seqId) = lr.get("ptnr1_label_asym_id", "ptnr1_label_seq_id");
			lr["ptnr1_label_atom_id"] = "ND2";
		}
		else
		{
			cif::tie(asymId, seqId) = lr.get("ptnr2_label_asym_id", "ptnr2_label_seq_id");
			lr["ptnr2_label_atom_id"] = "ND2";
		}
		
		for (auto r: db["atom_site"].find(
			cif::Key("label_asym_id") == asymId and
			cif::Key("label_comp_id") == "ASN" and
			cif::Key("label_seq_id") == seqId and
			(
				cif::Key("label_atom_id") == "OD1" or cif::Key("label_atom_id") == "ND2"
			)))
		{
			if (r["label_atom_id"].as<string>() == "ND2")
			{
				r["label_atom_id"] = "OD1";
				r["type_symbol"] = "O";

				auto a = db["atom_site_anisotrop"][cif::Key("id") == r["id"].as<string>()];
				if (not a.empty())
					a["type_symbol"] = "O";
			}
			else if (r["label_atom_id"].as<string>() == "OD1")
			{
				r["label_atom_id"] = "ND2";
				r["type_symbol"] = "N";

				auto a = db["atom_site_anisotrop"][cif::Key("id") == r["id"].as<string>()];
				if (not a.empty())
					a["type_symbol"] = "N";
			}
		}
		
		++numOfLinksFixed;
	}
	
	auto& nonO1 = dat["res_nonfix_O1"];
		
	// Fix O1-LINKs
	for (auto lr: db["struct_conn"].find(
		(cif::Key("ptnr1_label_atom_id") == "O1" and cif::Key("ptnr2_label_atom_id") == regex("C[2346]")) or
		(cif::Key("ptnr2_label_atom_id") == "O1" and cif::Key("ptnr1_label_atom_id") == regex("C[2346]"))))
	{
		string compId[2], atomId[2];
		cif::tie(compId[0], compId[1], atomId[0], atomId[1]) =
			lr.get("ptnr1_label_comp_id", "ptnr2_label_comp_id", "ptnr1_label_atom_id", "ptnr2_label_atom_id");
		
		if (not nonO1.find(cif::Key("id") == compId[0] or cif::Key("id") == compId[1]).empty())
			continue;
		
		if (VERBOSE)
			cerr << "Fixing non-standard carbohydrate LINK " << lr["id"] << endl;
		
		string asymId[2], seqId[2], entityId[2];
		
		smatch m;
		
		if (regex_match(atomId[0], m, regex("C([2346])")))
		{
			lr["ptnr1_label_atom_id"] = atomId[0] = "O" + m[1].str();
			lr["ptnr2_label_atom_id"] = "C1";
			
			cif::tie(asymId[0], seqId[0], asymId[1], seqId[1]) =
				lr.get("ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr2_label_asym_id", "ptnr2_label_seq_id");
		}
		else if (regex_match(atomId[1], m, regex("C([2346])")))
		{
			lr["ptnr1_label_atom_id"] = "C1";
			lr["ptnr2_label_atom_id"] = atomId[0] = "O" + m[1].str();
			
			cif::tie(asymId[0], seqId[0], asymId[1], seqId[1]) =
				lr.get("ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr1_label_asym_id", "ptnr1_label_seq_id");

			swap(compId[0], compId[1]);
		}

		// swap alt location
		auto alt = lr["ptnr1_label_alt_id"].as<string>();
		lr["ptnr1_label_alt_id"] = lr["ptnr2_label_alt_id"].as<string>();
		lr["ptnr2_label_alt_id"] = alt;
		
		// Update atom_site records
		for (auto a: db["atom_site"].find(
			cif::Key("label_asym_id") == asymId[1] and
			cif::Key("label_comp_id") == compId[1] and
			cif::Key("label_seq_id") == seqId[1] and
			cif::Key("label_atom_id") == "O1"))
		{
			a.setCascadeUpdate(false);

			a["label_atom_id"] = atomId[0];
			a["label_comp_id"] = compId[0];
			a["label_asym_id"] = asymId[0];
			a["label_entity_id"] = entityForAsym(asymId[0]);
			a["label_seq_id"] = seqId[0];
			
			// Also fix the pdb_auth fields (should this really be done? It is needed to export PDB again...)
			
			char chainID, iCode;
			int seqNum;
			
			tie(chainID, compId[0], seqNum, iCode) = MapBackResidue(db, asymId[0], compId[0], seqId[0]);
			
			a["auth_seq_id"] = to_string(seqNum);
			a["auth_comp_id"] = compId[0];
			a["auth_asym_id"] = string{ chainID };
			a["pdbx_PDB_ins_code"] = iCode == ' ' ? "" : string { iCode };
			a["auth_atom_id"] = atomId[0];
		}

		++numOfLinksFixed;
	}
	
	// Swap LINK labels to ensure directionality
	for (auto ll: db["struct_conn"].find(
		( cif::Key("ptnr1_label_atom_id") == "N" and cif::Key("ptnr2_label_atom_id") == "C") or
		( cif::Key("ptnr1_label_atom_id") == "P" and cif::Key("ptnr2_label_atom_id") == "O3'")))
	{
		for (string f: { "label_atom_id", "label_asym_id", "label_comp_id", "label_seq_id",
				"auth_asym_id", "auth_comp_id", "auth_seq_id", "symmetry" })
			SwapFields(ll, "ptnr1_" + f, "ptnr2_" + f);

		for (string f: { "label_alt_id", "PDB_ins_code" })
			SwapFields(ll, "pdbx_ptnr1_" + f, "pdbx_ptnr2_" + f);

		++numOfLinksFixed;
	}
	
	for (auto r: dat["link_swap"])
	{
		string atomId[2], compId[2];
		cif::tie(atomId[0], compId[0], atomId[1], compId[1]) =
			r.get("atom_id_1", "comp_id_1", "atom_id_2", "comp_id_2");
		
		for (auto ll: db["struct_conn"].find(
			cif::Key("ptnr1_label_atom_id") == atomId[0] and
			cif::Key("ptnr1_label_comp_id") == compId[0] and
			cif::Key("ptnr2_label_atom_id") == atomId[1] and
			cif::Key("ptnr2_label_comp_id") == compId[1]))
		{
			for (string f: { "label_atom_id", "label_asym_id", "label_comp_id", "label_seq_id",
					"auth_asym_id", "auth_comp_id", "auth_seq_id", "symmetry" })
				SwapFields(ll, "ptnr1_" + f, "ptnr2_" + f);
	
			for (string f: { "label_alt_id", "PDB_ins_code" })
				SwapFields(ll, "pdbx_ptnr1_" + f, "pdbx_ptnr2_" + f);

			++numOfLinksFixed;
		}
	}
	
	// // Set Refmac link types
	// for (auto t: dat["link_label"])
	// {
	// 	string atomId[2], compId[2], type;
	// 	cif::tie(atomId[0], compId[0], atomId[1], compId[1], type) =
	// 		t.get("atom_id_1", "comp_id_1", "atom_id_2", "comp_id_2", "refmac");
		
	// 	for (auto ll: db["struct_conn"].find(
	// 		cif::Key("ptnr1_label_atom_id") == atomId[0] and
	// 		cif::Key("ptnr1_label_comp_id") == compId[0] and
	// 		cif::Key("ptnr2_label_atom_id") == atomId[1] and
	// 		cif::Key("ptnr2_label_comp_id") == compId[1]))
	// 	{
	// 		ll["details"] = type;
	// 	}
	// }
	
	// change X atom_type to N for ASX/GLX 
	for (auto a: db["atom_site"].find(
		cif::Key("type_symbol") == "X" and
		(cif::Key("label_comp_id") == "ASX" or cif::Key("label_comp_id") == "GLX")))
	{
		a["type_symbol"] = "N";
		
		for (auto aa: db["atom_site_anisotrop"].find(cif::Key("id") == a["id"].as<string>()))
			aa["type_symbol"] = "N";
	}

	// Remove hydrogens and X
	for (auto h: db["atom_site"].find(cif::Key("type_symbol") == "H" or cif::Key("type_symbol") == "D" or cif::Key("type_symbol") == "X"))
	{
		if (VERBOSE)
			cerr << "Deleted (hetero) atom: " << h["id"] << endl;
		
		db["atom_site"].erase(h);

		++numOfAtomsDeleted;
	}
	
	// Remove unknown ligands when running in database mode.
	for (auto r: db["chem_comp"].find(cif::Key("id") == "UNL"))
	{
		// file contains UNL records
		if (serverMode)
		{
			if (c::CompoundFactory::instance().create("UNL") == nullptr)
				throw runtime_error("File contains unspecified UNL residues");

			if (VERBOSE)
				cerr << "File contains UNL records but also specifies the residue in the dictionary, accepting" << endl;

			break;
		}

		if (VERBOSE)
			cerr << "File contains UNL records but does not specify the residue in the dictionary, dropping" << endl;

		auto unls = db["atom_site"].find(cif::Key("label_comp_id") == "UNL");
		numOfAtomsDeleted += unls.size();
		
		for (auto a: unls)
			db["atom_site"].erase(a);

		db["pdbx_nonpoly_scheme"].erase(cif::Key("mon_id") == "UNL");
		db["pdbx_entity_nonpoly"].erase(cif::Key("comp_id") == "UNL");
		db["struct_conn"].erase(cif::Key("ptnr1_label_comp_id") == "UNL" or cif::Key("ptnr2_label_comp_id") == "UNL");
		db["chem_comp"].erase(r);

		break;
	}

//	if (serverMode)
//	{
//		for (auto a: db["atom_site"].find(cif::Key("label_comp_id") == "UNL"))
//		{
//			if (VERBOSE)
//				cerr << "Deleted (hetero) atom: " << a["id"] << endl;
//			
//			db["atom_site"].erase(cif::Key("label_comp_id") == "UNL");
//			
//			++numOfAtomsDeleted;
//		}
//	}
	
	for (auto a: db["atom_site"].find(cif::Key("label_comp_id") == "UNK"))
	{
		string atomId = a["label_atom_id"].as<string>();
		
		if (kBackBone.count(atomId) == 0 and atomId != "CB")
		{
			if (VERBOSE)
				cerr << "Deleted atom from unknown residue " << a["id"] << endl;
			
			db["atom_site"].erase(a);		// We have a crazy atom in an unknown residue

			++numOfAtomsDeleted;
		}
	}
	
	for (auto a: db["atom_site"].find(cif::Key("label_comp_id") == "GLY" and cif::Key("label_atom_id") == "CB"))
	{
		if (VERBOSE)
			cerr << "Deleted CB atom from GLY " << a["id"] << endl;
		
		db["atom_site"].erase(a);
		
		++numOfAtomsDeleted;
	}
	
	// Fix occupancy for Se in MSE when other side chain atoms have full occupancy
	for (auto a: db["atom_site"].find(
		cif::Key("label_comp_id") == "MSE" and
		cif::Key("label_alt_id") == "" and
		cif::Key("label_atom_id") == "SE" and
		cif::Key("occupancy") < 1.0))
	{
		string asymId, seqId;
		cif::tie(asymId, seqId) = a.get("label_asym_id", "label_seq_id");
		
		auto sideChainAtoms = db["atom_site"].find(
			cif::Key("label_asym_id") == asymId and
			cif::Key("label_seq_id") == seqId and
			cif::Key("label_comp_id") == "MSE" and
			(cif::Key("label_atom_id") == "CB" or cif::Key("label_atom_id") == "CG" or cif::Key("label_atom_id") == "CE") and
			cif::Key("occupancy") >= 1.0);
		
		if (not sideChainAtoms.empty())
		{
			if (VERBOSE)
				cerr << "Changed selenium occupancy to 1.00 for " << a["id"] << endl;
			
			a["occupancy"] = 1.0;
			
			++numOfOccupanciesReset;
		}
	}

	// delete aniso records when B factor < 2	
	for (auto a: db["atom_site"].find(cif::Key("B_iso_or_equiv") < 2.0))
	{
		string id = a["id"].as<string>();
		
		for (auto aa: db["atom_site_anisotrop"].find(cif::Key("id") == id))
		{
			db["atom_site_anisotrop"].erase(aa);
			
			if (VERBOSE)
				cerr << "Deleted aniso record for atom with id " << id << " since B-factor is less than 2" << endl;
			
			++numOfAnisosDeleted;
			
			a["B_iso_or_equiv"] = 2.0;
			++numOfBFactorsReset;
		}
	}
	
	// delete aniso records where one of the eigen values of U is less than 0
	set<string> delAniso;
	
	for (auto aa: db["atom_site_anisotrop"])
	{
		string id;
		float U[6];
		cif::tie(id, U[0], U[1], U[2], U[3], U[4], U[5]) =
			aa.get("id", "U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]");
		
		clipper::Matrix<> m(3, 3);
		m(0, 0) = U[0];
		m(1, 1) = U[1];
		m(2, 2) = U[2];
		m(0, 1) = m(1, 0) = U[3];
		m(0, 2) = m(2, 0) = U[4];
		m(1, 2) = m(2, 1) = U[5];
		
		if (m.eigen(true).front() < 0)
		{
			delAniso.insert(id);

			if (VERBOSE)
				cerr << "Deleted aniso record for atom with id " << id << " since one of the eigen values of U is less than zero" << endl;
			
			++numOfAnisosDeleted;
		}
	}

	for (auto id: delAniso)
		db["atom_site_anisotrop"].erase(cif::Key("id") == id);
	
	// Rename ATOM's with residue name WAT and atom id OW0 to HOH/O
	for (auto a: db["atom_site"].find(
		(cif::Key("label_comp_id") == "WAT" and cif::Key("label_atom_id") == "OW0") or
		(cif::Key("auth_comp_id") == "WAT" and cif::Key("auth_atom_id") == "OW0")))
	{
		a["label_comp_id"] = "HOH";
		a["label_atom_id"] = "O";
		a["auth_comp_id"] = "HOH";
		a["auth_atom_id"] = "O";

		++numOfRenamedWaters;
	}

	cout << endl
		 << "Report" << endl
		 << "Occupancies reset     : " << numOfOccupanciesReset << endl
		 << "B-factors reset       : " << numOfBFactorsReset << endl
		 << "Deleted (hetero) atoms: " << numOfAtomsDeleted << endl
		 << "Deleted ANISOU records: " << numOfAnisosDeleted << endl
		 << "Deleted LINK   records: " << numOfLinksDeleted << endl
		 << "Added   LINK   records: " << numOfLinksAdded << endl
		 << "Fixed   LINK   records: " << numOfLinksFixed << endl
		 << "Renamed water atoms   : " << numOfRenamedWaters << endl
		 ;
		 
	if (vm.count("output"))
	{
		db.add_software("prepper", "other", get_version_nr(), get_version_date());

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
		
		pdb.save(vm["output"].as<string>());
	}

	return result;
}
