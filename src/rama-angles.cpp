/* 
   Created by: Maarten L. Hekkelman
   Date: woensdag 09 januari, 2019
*/

#include "pdb-redo.h"

#include "cif++/Structure.h"
#include "cif++/Secondary.h"

#include <zeep/http/webapp.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>

#include "ramachandran.h"

#include "mrsrc.h"

#include "pr-server.hpp"

using namespace std;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace zh = zeep::http;
namespace ba = boost::algorithm;

#if LIBZEEP_VERSION_MAJOR >= 4
#include <zeep/el/element.hpp>

namespace el = zeep::el;
using json = el::element;

#else
#include <nlohmann/json.hpp>
#include <zeep/http/webapp/el.hpp>
namespace el = zeep::http::el;

using json = nlohmann::json;
#endif

using mmcif::Structure;

// --------------------------------------------------------------------

fs::path
	gPDBDIR = "/srv/data/pdb/mmCIF/",
	gPDBREDODIR = "/srv/pdb_redo/";

fs::path gExePath;

const char* PID_FILE = "/var/run/rama-angles.pid";

// --------------------------------------------------------------------

struct ResidueInfo
{
	ResidueInfo(const ResidueInfo&) = default;
	ResidueInfo& operator=(const ResidueInfo&) = default;

	ResidueInfo(const string& residueName, int residueNumber, int authorResidueNumber, const string& authorInsertionCode)
		: residue_name(residueName), residue_number(residueNumber)
		, author_residue_number(authorResidueNumber), author_insertion_code(authorInsertionCode) {}

	ResidueInfo(const string& residueName, int residueNumber, int authorResidueNumber, const string& authorInsertionCode, float psi, float phi, bool cis = false)
		: residue_name(residueName), residue_number(residueNumber)
		, author_residue_number(authorResidueNumber), author_insertion_code(authorInsertionCode)
		, psi(psi), phi(phi), cis(cis) {}
	
	string	residue_name;
	int		residue_number;
	int		author_residue_number;
	string	author_insertion_code;
	float	psi = 360, phi = 360;
	string	rama;
	
//	string	rota;
//	string	cis_peptide;

	bool	cis = false;

	json toJSON() const;
};

json ResidueInfo::toJSON() const
{
	json result = {
		{ "psi", nullptr },
//		{ "cis_peptide", nullptr },
		{ "residue_number", residue_number },
		{ "author_residue_number", author_residue_number },
		{ "rama", nullptr },
		{ "phi", nullptr },
		{ "author_insertion_code", author_insertion_code },
		{ "residue_name", residue_name },
//		{ "rota", rota }
	};
	
	if (phi < 360 and psi < 360)
	{
		result["phi"] = phi;
		result["psi"] = psi;
		result["rama"] = rama;
	}
	
	if (cis)
		result["cis_peptide"] = "Y";
	
	return result;
}

// --------------------------------------------------------------------

struct ChainInfo
{
	ChainInfo(const ChainInfo&) = default;
	ChainInfo& operator=(const ChainInfo&) = default;
	
	ChainInfo(const string& entityID, int modelID, const string& chainID, const string& asymID)
		: entity_id(entityID), model_id(modelID), chain_id(chainID), struct_asym_id(asymID) {}
	
	string				entity_id;
	int					model_id;
	string				chain_id;
	string				struct_asym_id;
	vector<ResidueInfo>	residues;
};

vector<ChainInfo> CalculateChainInfo(Structure& structure)
{
	using namespace mmcif;
	
	vector<ChainInfo> result;

	for (auto& poly: structure.polymers())
	{
		string asymID = poly.asymID();
		
		result.emplace_back(poly.entityID(), 1, poly.chainID(), asymID);
		auto& residues = result.back().residues;
		
		for (auto& r: poly)
			residues.emplace_back(r.compoundID(), r.seqID(), stoi(r.authSeqID()), r.authInsCode(), r.psi(), r.phi(), r.isCis());
		
		for (size_t i = 1; i + 1 < residues.size(); ++i)
		{
			auto& ri = residues[i];

			bool prePro = (i + 1 < residues.size() and residues[i + 1].residue_name == "PRO" and residues[i].residue_name != "PRO");
			
			switch (calculateRamachandranScore(ri.residue_name, prePro, ri.phi, ri.psi))
			{
				case rsNotAllowed:	ri.rama = "OUTLIER"; break;
				case rsAllowed:		ri.rama = "Allowed"; break;
				case rsFavoured:	ri.rama = "Favored"; break;
			}
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

json CreateJSONForStructureFile(const string& inFile)
{
	if (VERBOSE)
		cerr << "Reading file " << inFile << endl;

	mmcif::File f(inFile);

	if (VERBOSE)
		cerr << "Creating structure" << endl;

	Structure structure(f);
	string id = f.data().getName();
	ba::to_lower(id);

	// --------------------------------------------------------------------

	auto ciList = CalculateChainInfo(structure);
	set<string> entities;
	for (auto& ci: ciList)
		entities.insert(ci.entity_id);

	// --------------------------------------------------------------------

	json molecules;
	
	for (string entityID: entities)
	{
		json chains;
		
		for (auto& ci: ciList)
		{
			if (ci.entity_id != entityID)
				continue;
			
			json residues;
			
			for (auto& ri: ci.residues)
				residues.emplace_back(ri.toJSON());
			
			json chain = {
				{
					"models", {
						{
							{ "model_id", 1 },
							{ "residues", move(residues) }
						}
					}
				},
				{ "chain_id", ci.chain_id },
				{ "struct_asym_id", ci.struct_asym_id }
			};
			
			chains.emplace_back(move(chain));
		}
		
		json molecule = {
			{ "entity_id", entityID },
			{ "chains", move(chains) }
		};
		
		molecules.push_back(move(molecule));
	}

	json result =
	{
		{
			id, {
				{ "molecules", move(molecules) }
			}
		}
	};
	
	return result;
}

// --------------------------------------------------------------------

const char kRamaAnglesNS[] = "http://pdb-redo.eu/ns/rama-angles";

class RamaAnglesServer : public zh::webapp
{
  public:
	RamaAnglesServer();
	~RamaAnglesServer();
	
    void handle_rama_request(const zh::request& request, const el::scope& scope, zh::reply& reply);
	virtual void load_template(const string& file, zeep::xml::document& doc);
	void handle_file(const zh::request& request, const el::scope& scope, zh::reply& reply);
};

RamaAnglesServer::RamaAnglesServer()
	: webapp(kRamaAnglesNS)
{
	mount("orig",		std::bind(&RamaAnglesServer::handle_rama_request, this, placeholders::_1, placeholders::_2, placeholders::_3));
	mount("redo",		std::bind(&RamaAnglesServer::handle_rama_request, this, placeholders::_1, placeholders::_2, placeholders::_3));
	mount("style.css",	std::bind(&RamaAnglesServer::handle_file, this, placeholders::_1, placeholders::_2, placeholders::_3));
}

RamaAnglesServer::~RamaAnglesServer()
{
}

void RamaAnglesServer::handle_rama_request(const zh::request& request, const el::scope& scope, zh::reply& reply)
{
	fs::path file = scope["baseuri"].as<string>();

	string selector = file.parent_path().string();
	string pdbid = file.filename().string();
	
	ba::to_lower(pdbid);
	
	if (pdbid.length() != 4)	// sja... TODO!!! fix
		throw runtime_error("not a valid pdbid");
	
	fs::path f;
	if (selector == "orig")
		f = gPDBDIR / pdbid.substr(1, 2) / (pdbid + ".cif.gz");
	else if (selector == "redo")
		f = gPDBREDODIR / pdbid.substr(1, 2) / pdbid / (pdbid + "_final.cif");
	
	if (not fs::exists(f))
		throw runtime_error("pdb file not found: " + f.string());

	json result = CreateJSONForStructureFile(f.string());
	
#if LIBZEEP_VERSION_MAJOR >= 4
	reply.set_content(result);
#else
	reply.set_content(result.dump(), "application/json");
#endif
}

void RamaAnglesServer::load_template(const std::string& file, zeep::xml::document& doc)
{
	mrsrc::rsrc rsrc(file);
	if (not rsrc)
		throw runtime_error("missing template");
	
	string data(rsrc.data(), rsrc.size());
	doc.read(data);
}

void RamaAnglesServer::handle_file(const zh::request& request, const el::scope& scope, zh::reply& reply)
{
	using namespace boost::local_time;
	using namespace boost::posix_time;

	fs::path file = scope["baseuri"].as<string>();

	mrsrc::rsrc rsrc(file.string());

	if (not rsrc)
		create_error_reply(request, zh::not_found, "The requested file was not found on this 'server'", reply);
	else
	{
		// compare with the date/time of our executable, since we're reading resources :-)
		string ifModifiedSince;
		for (const zeep::http::header& h : request.headers)
		{
			if (ba::iequals(h.name, "If-Modified-Since"))
			{
				local_date_time modifiedSince(local_sec_clock::local_time(time_zone_ptr()));

				local_time_input_facet* lif1(new local_time_input_facet("%a, %d %b %Y %H:%M:%S GMT"));

				stringstream ss;
				ss.imbue(std::locale(std::locale::classic(), lif1));
				ss.str(h.value);
				ss >> modifiedSince;

				local_date_time fileDate(from_time_t(fs::last_write_time(gExePath)), time_zone_ptr());

				if (fileDate <= modifiedSince)
				{
					reply = zeep::http::reply::stock_reply(zeep::http::not_modified);
					return;
				}

				break;
			}
		}

		string data(rsrc.data(), rsrc.size());
		string mimetype = "text/plain";
	
		if (file.extension() == ".css")
			mimetype = "text/css";
		else if (file.extension() == ".js")
			mimetype = "text/javascript";
		else if (file.extension() == ".png")
			mimetype = "image/png";
		else if (file.extension() == ".svg")
			mimetype = "image/svg+xml";
		else if (file.extension() == ".html" or file.extension() == ".htm")
			mimetype = "text/html";
		else if (file.extension() == ".xml" or file.extension() == ".xsl" or file.extension() == ".xslt")
			mimetype = "text/xml";
		else if (file.extension() == ".xhtml")
			mimetype = "application/xhtml+xml";
	
		reply.set_content(data, mimetype);
	
		local_date_time t(local_sec_clock::local_time(time_zone_ptr()));
		local_time_facet* lf(new local_time_facet("%a, %d %b %Y %H:%M:%S GMT"));
		
		stringstream s;
		s.imbue(std::locale(std::cout.getloc(), lf));
		
		ptime pt = from_time_t(boost::filesystem::last_write_time(gExePath));
		local_date_time t2(pt, time_zone_ptr());
		s << t2;
	
		reply.set_header("Last-Modified", s.str());
	}
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " [coordinatesfile [outputfile]] options");
	visible_options.add_options()
		("xyzin",				po::value<string>(),	"coordinates file")
		("output,o",			po::value<string>(),	"Write output to this file instead of stdout")
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		
		("server",										"Start as web service")
		("address",				po::value<string>(),	"External address, default is 0.0.0.0")
		("port",				po::value<uint16_t>(),	"Port to listen to, default is 10336")
		("no-daemon,F",									"Do not fork into background")
		("user,u",				po::value<string>(),	"User to run the daemon")
		("logfile",				po::value<string>(),	"Logfile to write to, default /var/log/rama-angles.log")

		("pdb-dir",				po::value<string>(),	"PDB mmCIF data directory")
		("pdb-redo-dir",		po::value<string>(),	"PDB-REDO data directory")
		("ccp4-dir",			po::value<string>(),	"Directory where ccp4 is installed")

		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		cerr << visible_options << endl;
		exit(0);
	}
	
	if (vm.count("xyzin") == vm.count("server"))
	{
		cerr << visible_options << endl;
		exit(1);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();

	if (vm.count("pdb-dir"))
		gPDBDIR = vm["pdb-dir"].as<string>();

	if (vm.count("pdb-redo-dir"))
		gPDBREDODIR = vm["pdb-redo-dir"].as<string>();

	const char* ccp4 = getenv("CCP4");
	string ccp4Dir = ccp4 ? ccp4 : "";
	if (vm.count("ccp4-dir"))
		ccp4Dir = vm["ccp4-dir"].as<string>();
	if (not ccp4Dir.empty())
	{
		auto clibdMon = fs::path(ccp4Dir) / "lib" / "data" / "monomers";
		setenv("CLIBD_MON", clibdMon.string().c_str(), 1);
	}

	// --------------------------------------------------------------------

	if (vm.count("server"))
	{
		char exePath[PATH_MAX + 1];
		int r = readlink("/proc/self/exe", exePath, PATH_MAX);
		if (r > 0)
		{
			exePath[r] = 0;
			gExePath = fs::system_complete(exePath);
		}
		
		if (not fs::exists(gExePath))
			gExePath = fs::system_complete(argv[0]);
	
		string address = "0.0.0.0";
		if (vm.count("address"))
			address = vm["address"].as<string>();
		
		uint16_t port = 10336;
		if (vm.count("port"))
			port = vm["port"].as<uint16_t>();

		string user = "nobody";
		if (vm.count("user"))
			user = vm["user"].as<string>();
		
		bool fork = vm.count("no-daemon") == 0;

		if (fork)
			Daemonize();
		
		string logFile = "/var/log/rama-angles.log";
		if (vm.count("logfile"))
			logFile = vm["logfile"].as<string>();

		auto serverFactory = []() -> zeep::http::server*
		{
			return new RamaAnglesServer();
		};
		
		for (;;)
		{
			if (fork and not logFile.empty())
				OpenLogFile(logFile);
			
			if (RunMainLoop(address, port, user, serverFactory))
				continue;
			
			break;
		}
		
		return 0;
	}

	// --------------------------------------------------------------------

	json result = CreateJSONForStructureFile(vm["xyzin"].as<string>());

	ofstream of;
	if (vm.count("output"))
	{
		of.open(vm["output"].as<string>());
		if (not of.is_open())
		{
			cerr << "Could not open output file" << endl;
			exit(1);
		}
		
		of << result << endl;
	}
	else
		cout << setw(2) << result << endl;

	return 0;	
}
