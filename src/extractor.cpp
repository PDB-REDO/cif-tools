#include "pdb-redo.h"

#include <sys/wait.h>

#include <fstream>
#include <chrono>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>


#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include "cif++/Cif++.h"
#include "cif++/PDB2Cif.h"
#include "cif++/Cif2PDB.h"
#include "cif++/CifUtils.h"
#include "cif++/Structure.h"
#include "cif++/TlsParser.h"
#include "cif++/ResolutionCalculator.h"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace c = mmcif;

using cif::iequals;

// --------------------------------------------------------------------

struct s1
{
	float val;
	int w, p;
	
	s1(cif::Row& db, const char* fld, int w = 9, int p = 4)
		: val(db[fld].as<float>()), w(w), p(p) {}
};

ostream& operator<<(ostream& os, const s1& v)
{
	os << fixed << right << setw(v.w) << setprecision(v.p) << v.val;
	return os;
};

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	int result = 0;
	
	po::options_description visible_options("extractor " + VERSION_STRING + " options coordinatesfile structurefile" );
	visible_options.add_options()
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("output,o",	po::value<string>(),	"The output file, default is stdout")
		("tls-out",		po::value<string>(),	"Name for the TLS file to create, default is calculated based on input file")
		("pdb-redo-data",po::value<string>(),	"The PDB-REDO dat file" /*, default is the built in one"*/)
		("dict",		po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
		;

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("coordinates",	po::value<string>(),	"Coordinates file")
		("reflections",	po::value<string>(),	"Reflections file")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("coordinates", 1);
	p.add("reflections", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("coordinates") == 0 or vm.count("reflections") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}
	
	fs::path pdbRedoDataFile = "pdb-redo-data.cif";
	if (vm.count("pdb-redo-data"))
		pdbRedoDataFile = vm["pdb-redo-data"].as<string>();
	std::ifstream pdbRedoData(pdbRedoDataFile);
	if (not pdbRedoData.is_open())
		throw runtime_error("Could not open pdb-redo-data file");
	cif::File pdbRedoDataCif(pdbRedoData);
	auto& redoDB = pdbRedoDataCif["PDB_REDO_DAT"];

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	// Load dict, if any
	
	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());

	// Input coordinates file
	
	fs::path coordinates = vm["coordinates"].as<string>();
	c::File coordFile(coordinates);
	c::Structure structure(coordFile, 1);
	auto& coord = coordFile.data();
	
	string entryId = coord["entry"].front()["id"].as<string>();
	if (entryId.empty())
	{
		coordFile.save("/tmp/dump.cif");
		throw runtime_error("Missing _entry.id in coordinates file");
	}
	
	// Input reflections file
	
	fs::path reflections = vm["reflections"].as<string>();
	cif::File reflnFile;
	reflnFile.load(reflections);
	auto& refln = reflnFile.firstDatablock();
	
	// Output file, create it
	fs::path coordinatesFileName = coordinates.filename();
	if (coordinatesFileName.extension() == ".gz" or coordinatesFileName.extension() == ".bz2")
		coordinatesFileName = coordinatesFileName.stem();
	
	fs::path outputStem = coordinatesFileName.stem();
	if (vm.count("output") > 0)
		outputStem = vm["output"].as<string>();
	
	fs::path output = outputStem;
	output += ".extracted";
	fs::ofstream outputFile(output);
	
	fs::path tlsout = outputStem;
	tlsout += ".tls";
	
	if (vm.count("tls-out"))
		tlsout = vm["tls-out"].as<string>();
	
	fs::ofstream tlsoutFile(tlsout);
	
	cif::File extractedFile;
	extractedFile.append(new cif::Datablock(entryId));
	auto& extractedDb = extractedFile.firstDatablock();
	
	// Ready to do some work
	
	double a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = coord["cell"][cif::Key("entry_id") == entryId]
		.get("length_a", "length_b", "length_c",
			 "angle_alpha", "angle_beta", "angle_gamma");
	
	extractedDb["cell"].emplace({
			{ "length_a",    a },
			{ "length_b",    b },
			{ "length_c",    c },
			{ "angle_alpha", alpha },
			{ "angle_beta",  beta },
			{ "angle_gamma", gamma }
		});

	auto refine = coord["refine"][cif::Key("entry_id") == entryId];
	if (refine.empty())
		throw runtime_error("No refinement data found");
	
	double reso, hires = 99, lowres = 0;
	cif::tie(reso, hires, lowres) = refine.get("ls_d_res_high", "ls_d_res_high", "ls_d_res_low");
	
	string spacegroup = coord["symmetry"]
		[cif::Key("entry_id") == entryId]
		["space_group_name_H-M"].as<string>();
	
	if (spacegroup == "P 1-")
		spacegroup = "P -1";
	else if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw runtime_error("No spacegroup, cannot continue");
	
	float wavelength = 0;
	for (auto s: coord["diffrn_source"])
	{
		wavelength = s["wavelength"].as<float>();
		break;
	}

	if (wavelength == 0)
	{
		for (auto s: coord["diffrn_radiation_wavelength"])
		{
			wavelength = s["wavelength"].as<float>();
			break;
		}
	}
	
	// Read reflections file to calculate resolution low and high
	c::ResolutionCalculator rc(a, b, c, alpha, beta, gamma);
	uint32_t refcnt = 0, tstcnt = 0;
	
	for (auto r: refln["refln"])
	{
		int h, k, l;
		char flag;
		
		cif::tie(h, k, l, flag) = r.get("index_h", "index_k", "index_l", "status");
		
		double res = rc(h, k, l);
		
		if (hires > res)
			hires = res;
		if (lowres < res)
			lowres = res;
		
		++refcnt;
		if (flag == 'f')
			++tstcnt;
	}
	
	// Find R-factor and R-free
	float rfactor = 0.999, rfree = 0.999;

	for (auto fld: { "_refine.ls_R_factor_R_work", "_refine.ls_R_factor_obs",
		"_pdbx_refine.R_factor_obs_no_cutoff", "_pdbx_refine.R_factor_obs_4sig_cutoff" })
	{
		string category, tag;
		tie(category, tag) = cif::splitTagName(fld);
		
		auto r = coord[category][cif::Key("entry_id") == entryId];
		
		if (r[tag].empty())
			continue;
		
		rfactor = r[tag].as<float>();
		break;
	}
		
	for (auto fld: { "_refine.ls_R_factor_R_free" })
	{
		string category, tag;
		tie(category, tag) = cif::splitTagName(fld);
		
		auto r = coord[category][cif::Key("entry_id") == entryId];
		
		if (r[tag].empty())
			continue;
		
		rfree = r[tag].as<float>();
		break;
	}
		
	// Calculate the average B-factor
	double bTotal = 0;
	uint32_t count = 0;
	for (auto r: coord["atom_site"])
	{
		if (not r["B_iso_or_equiv"].empty())
		{
			++count;
			bTotal += r["B_iso_or_equiv"].as<double>();
		}
	}
	
	double bAvg = bTotal / count;
	
	// refinement program
	string program;
	for (auto p: coord["software"].find(cif::Key("classification") == "refinement"))
	{
		program = p["name"].as<string>();
		break;
	}

	// Store other info
	extractedDb["extracted_info"].emplace({
		{ "resolution", reso },
		{ "resolution_low", lowres },
		{ "resolution_high", hires },
		{ "spacegroup", spacegroup },
		{ "wavelength", wavelength },
		{ "R-factor", rfactor },
		{ "R-free", rfree },
		{ "B-avg", bAvg },
		{ "reflection-count", refcnt },
		{ "testset-count", tstcnt },
		{ "testset-perc", (100.0 * tstcnt) / refcnt },
		{ "refinement-program", program },
		{ "protein", coord["entity_poly"].exists(cif::Key("type") == "polypeptide(L)" or cif::Key("type") == "polypeptide(D)") },
		{ "dna", coord["entity_poly"].exists(cif::Key("type") == "polydeoxyribonucleotide") },
		{ "rna", coord["entity_poly"].exists(cif::Key("type") == "polyribonucleotide") },
		{ "polysaccharide", coord["entity_poly"].exists(cif::Key("type") == "polysaccharide(L)" or cif::Key("type") == "polysaccharide(D)") },
	});

	// Parse TLS group specifications. If needed...
	
	for (auto tls: coord["pdbx_refine_tls"])
	{
		string refineTlsId = tls["id"].as<string>();
		
		stringstream tlsGroup;

		bool hasRanges = false;
		tlsGroup << "TLS    Group " << setw(3) << right << refineTlsId << endl; 

		for (auto t: coord["pdbx_refine_tls_group"].find(cif::Key("refine_tls_id") == refineTlsId))
		{
			string selectionDetails = t["selection_details"].as<string>();

			if (not selectionDetails.empty())
			{
				auto s = cif::ParseSelectionDetails(program, selectionDetails);
				
				if (s == nullptr)	// give up if we could not parse this.
					continue;
				
				for (auto r: s->GetRanges(coord, true))
				{
					string chainID;
					int from, to;
					tie(chainID, from, to) = r;
					
					tlsGroup << "RANGE  '" << chainID << setw(4) << from << ".' '" << chainID << setw(4) << to << ".' ALL" << endl;

					hasRanges = true;
				}
				
				continue;
			}
			
			// collect the auth based range, not the label one. 
			
			string begAsymId, endAsymId;
			int begSeqId, endSeqId;
			
			cif::tie(begAsymId, begSeqId, endAsymId, endSeqId)
				= t.get("beg_auth_asym_id", "beg_auth_seq_id", "end_auth_asym_id", "end_auth_seq_id");
				
			tlsGroup << "RANGE  '" << begAsymId << setw(4) << begSeqId << ".' '" << endAsymId << setw(4) << endSeqId << ".' ALL" << endl;
			hasRanges = true;
		}
		
		if (not hasRanges)
		{
			if (cif::VERBOSE)
				cout << "Skipping TLS group " << refineTlsId << " since selection is empty" << endl;
			
			continue;
		}
		
		float s22ms11 = tls["S[2][2]"].as<float>() - tls["S[1][1]"].as<float>();
		float s11ms33 = tls["S[1][1]"].as<float>() - tls["S[3][3]"].as<float>();
		
		tlsoutFile << tlsGroup.str();
		
		tlsoutFile << "ORIGIN " << s1(tls, "origin_x") << ' ' << s1(tls, "origin_y") << ' ' << s1(tls, "origin_z") << endl
				   << "T   " << s1(tls, "T[1][1]") << ' ' << s1(tls, "T[2][2]") << ' ' << s1(tls, "T[3][3]") << ' ' << s1(tls, "T[1][2]") << ' ' << s1(tls, "T[1][3]") << ' ' << s1(tls, "T[2][3]") << endl
				   << "L   " << s1(tls, "L[1][1]") << ' ' << s1(tls, "L[2][2]") << ' ' << s1(tls, "L[3][3]") << ' ' << s1(tls, "L[1][2]") << ' ' << s1(tls, "L[1][3]") << ' ' << s1(tls, "L[2][3]") << endl 
				   << "S   " << fixed << right << setw(9) << setprecision(4) << s22ms11 << ' '
				   			 << fixed << right << setw(9) << setprecision(4) << s11ms33 << ' '
				   			 << s1(tls, "S[1][2]") << ' ' << s1(tls, "S[1][3]") << ' ' << s1(tls, "S[2][3]") << ' '
				   			 << s1(tls, "S[2][1]") << ' ' << s1(tls, "S[3][1]") << ' ' << s1(tls, "S[3][2]")
				   			 << endl
				   << endl;
	}
	
	const cif::iset kBackBoneAtoms = { "N", "CA", "C", "O", "OXT" };
	
	// No build list
	for (auto l: coord["struct_conn"])
	{
		for (string prefix: { "ptnr1_", "ptnr2_" })
		{
			string atomId, asymId, monId, authSeqId;
			int seqId;

			cif::tie(atomId, asymId, monId, seqId, authSeqId) = l.get(
				(prefix + "label_atom_id").c_str(),
				(prefix + "label_asym_id").c_str(),
				(prefix + "label_comp_id").c_str(),
				(prefix + "label_seq_id").c_str(),
				(prefix + "auth_seq_id").c_str());
	
			if (not c::kAAMap.count(monId) and not iequals(monId, "HOH"))
				continue;

			string atomId2 = l[prefix == "ptnr1_" ? 
					"ptnr2_label_atom_id" : "ptnr1_label_atom_id"].as<string>();
	
			string special;
			
			if (iequals(monId, "HOH"))
				special = "H2O-keep";
			else if (kBackBoneAtoms.count(atomId) == 0)
				special = "no-build";
			else if (iequals(atomId, "N") and not iequals(atomId, "C"))
				special = "bbn-keep";
			else if (not iequals(atomId, "C"))
				special = "bbo-keep"; 
			
			if (special.empty())
				continue;
			
			string pdbStrandId, pdbMonId, pdbInsCode;
			int pdbSeqNum;

			tie(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode) =
				structure.MapLabelToPDB(asymId, seqId, monId, authSeqId);
			
			extractedDb["res_info"].emplace({
				{ "type", special },

				{ "asym_id", asymId },
				{ "mon_id", monId },
				{ "seq_id", seqId },
				
				{ "pdb_strand_id", pdbStrandId },
				{ "pdb_seq_num", pdbSeqNum },
				{ "pdb_mon_id", pdbMonId },
				{ "pdb_ins_code", pdbInsCode }
			});
		}		
	}
	
	vector<cif::Row> residues;
	cif::iset resIgnore;
	for (auto r: redoDB["res_hetero_ignore"])
		resIgnore.insert(r["name"].as<string>());
	
	// Write out ligand list
	for (auto r: coord["pdbx_poly_seq_scheme"])
		residues.push_back(r);

	for (auto r: coord["pdbx_nonpoly_scheme"])
		residues.push_back(r);

	for (auto r: residues)
	{
		string asymId, seqId, monId, pdbStrandId, pdbMonId, pdbInsCode;
		int pdbSeqNum;

		cif::tie(asymId, monId, seqId, pdbStrandId, pdbMonId, pdbSeqNum, pdbInsCode) =
			r.get("asym_id", "mon_id", "seq_id", "pdb_strand_id", "pdb_mon_id", "pdb_seq_num", "pdb_ins_code");

		if (c::CompoundFactory::instance().isKnownPeptide(monId) or resIgnore.count(monId))
			continue;
		
		extractedDb["ligand"].emplace({
			{ "asym_id", asymId },
			{ "mon_id", monId },
			{ "seq_id", seqId },
			
			{ "pdb_strand_id", pdbStrandId },
			{ "pdb_seq_num", pdbSeqNum },
			{ "pdb_mon_id", pdbMonId },
			{ "pdb_ins_code", pdbInsCode }
		});
	}
						
	// Write out results
	extractedFile.save(outputFile);
	
	return result;
}

