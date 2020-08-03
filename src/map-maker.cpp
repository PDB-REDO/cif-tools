// my humble attempt to create map files

#include "pdb-redo.h"

#include <clipper/clipper-ccp4.h>

#include <zeep/el/element.hpp>
#include <zeep/http/webapp.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>

#include "cif++/MapMaker.h"

#include "cif++/mrsrc.h"
#include "pr-server.hpp"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace c = mmcif;
namespace zh = zeep::http;
namespace el = zeep::el;
namespace ba = boost::algorithm;

// --------------------------------------------------------------------

fs::path
	gPDBDIR = "/srv/data/pdb/mmCIF/",
	gPDBREDODIR = "/srv/pdb_redo/";

fs::path gExePath;

const char* PID_FILE = "/var/run/map-maker.pid";

// --------------------------------------------------------------------

class self_destructing_file : public fs::ifstream
{
  public:
	self_destructing_file(fs::path p)
		: fs::ifstream(p, ios::binary), m_path(p) {}
	
	virtual ~self_destructing_file()
	{
		if (fs::exists(m_path))
			fs::remove(m_path);
	}

  private:
	fs::path m_path;
};

// --------------------------------------------------------------------

ostream& operator<<(ostream& os, const clipper::Coord_grid& c)
{
	os << '{' << c[0] << ',' << c[1] << ',' << c[2] << '}';
	return os;
}

// --------------------------------------------------------------------

void carveOutNXMap(const mmcif::Structure& structure, clipper::Xmap<float>& xMap, clipper::NXmap<float> nxMap)
{
	const float kMargin = 5;

	// bepaal eerst de box om het model met een marge van kMargin
	clipper::Coord_orth oMin, oMax;
	oMin[0] = oMin[1] = oMin[2] = numeric_limits<int>::max();
	oMax[0] = oMax[1] = oMax[2] = numeric_limits<int>::min();
	
	for (auto& atom: structure.atoms())
	{
		auto l = atom.location();
		if (oMin[0] > l.mX)
			oMin[0] = l.mX;
		if (oMax[0] < l.mX)
			oMax[0] = l.mX;
		if (oMin[1] > l.mY)
			oMin[1] = l.mY;
		if (oMax[1] < l.mY)
			oMax[1] = l.mY;
		if (oMin[2] > l.mZ)
			oMin[2] = l.mZ;
		if (oMax[2] < l.mZ)
			oMax[2] = l.mZ;
	}
	
	oMin[0] -= kMargin;
	oMax[0] += kMargin;
	oMin[1] -= kMargin;
	oMax[1] += kMargin;
	oMin[2] -= kMargin;
	oMax[2] += kMargin;

	clipper::Coord_frac fMin{oMin.coord_frac(xMap.cell())}, fMax{oMax.coord_frac(xMap.cell())};
	clipper::Coord_grid gMin{fMin.coord_grid(xMap.grid_sampling())}, gMax{fMax.coord_grid(xMap.grid_sampling())};

	// Definieer nu de nieuwe nxMap
	// nxMap.init(xMap.cell(), xMap.grid_sampling(), clipper::Grid_range{gMin, gMax});

	// construct an nxmap with a 

	// nxMap.init()


	// // copy data
	// auto i0 = clipper::Xmap_base::Map_reference_coord(xMap, gMin);
	// for (auto iu = i0; iu.coord().u() <= gMax[0]; iu.next_u())
	// 	for (auto iv = iu; iv.coord().v() <= gMax[1]; iv.next_v())
	// 		for (auto iw = iv; iw.coord().w() <= gMax[2]; iw.next_w())
	// 		{
	// 			auto c = iw.coord();
	// 			nxMap.set_data(c - gMin, xMap.get_data(c));
	// 		}
}

// --------------------------------------------------------------------

clipper::Grid_range calculateRange(const mmcif::Structure& structure, clipper::Xmap<float>& xMap, float border)
{
	// bepaal eerst de box om het model met een marge van kMargin
	clipper::Coord_orth oMin, oMax;
	oMin[0] = oMin[1] = oMin[2] = numeric_limits<int>::max();
	oMax[0] = oMax[1] = oMax[2] = numeric_limits<int>::min();
	
	for (auto& atom: structure.atoms())
	{
		auto l = atom.location();
		if (oMin[0] > l.mX)
			oMin[0] = l.mX;
		if (oMax[0] < l.mX)
			oMax[0] = l.mX;
		if (oMin[1] > l.mY)
			oMin[1] = l.mY;
		if (oMax[1] < l.mY)
			oMax[1] = l.mY;
		if (oMin[2] > l.mZ)
			oMin[2] = l.mZ;
		if (oMax[2] < l.mZ)
			oMax[2] = l.mZ;
	}
	
	oMin[0] -= border;
	oMax[0] += border;
	oMin[1] -= border;
	oMax[1] += border;
	oMin[2] -= border;
	oMax[2] += border;

	clipper::Coord_frac fMin{oMin.coord_frac(xMap.cell())}, fMax{oMax.coord_frac(xMap.cell())};
	clipper::Coord_grid gMin{fMin.coord_grid(xMap.grid_sampling())}, gMax{fMax.coord_grid(xMap.grid_sampling())};

	return { gMin, gMax };
}

// --------------------------------------------------------------------

void createMap(fs::path xyzin, fs::path hklin, string type,
	bool masked, float border, float samplingRate, zh::reply& reply)
{
	fs::path mapout = fs::unique_path("map-maker-%%%%-%%%%-%%%%-%%%%.map");
	
	mmcif::File f(xyzin);
	mmcif::Structure structure(f);

	mmcif::MapMaker<float> mm;
	
	mm.loadMTZ(hklin, samplingRate);

	// what range?
	auto range = masked ? calculateRange(structure, mm.fb(), border) : mm.fb().get().grid_asu();

	if (type == "density")
		mm.fb().write_masked(mapout, range);
	else if (type == "difference")
		mm.fd().write_masked(mapout, range);
	else if (type == "anomalous")
	{
		if (mm.fa().get().is_null())
			throw runtime_error("Anomalous map not present");
		mm.fa().write_masked(mapout, range);
	}
	else
		throw runtime_error("Unsupported map type " + type);
	
	reply.set_content(new self_destructing_file(mapout), "application/octet-stream");
}

// --------------------------------------------------------------------

const char kMapMakerNS[] = "http://pdb-redo.eu/ns/map-maker";

class MapMakerServer : public zh::webapp
{
  public:
	MapMakerServer();
	~MapMakerServer();
	
    void handle_map_request(const zh::request& request, const el::scope& scope, zh::reply& reply);
	virtual void load_template(const string& file, zeep::xml::document& doc);
	void handle_file(const zh::request& request, const el::scope& scope, zh::reply& reply);
};

MapMakerServer::MapMakerServer()
	: webapp(kMapMakerNS)
{
	// create our temp directory and move there
	fs::path tmpDir = fs::temp_directory_path() / ("map-maker-" + to_string(getpid()));
	fs::create_directories(tmpDir);
	fs::current_path(tmpDir);

	using namespace std::placeholders;

	mount("map",		&MapMakerServer::handle_map_request);
	mount("style.css",	&MapMakerServer::handle_file);
}

MapMakerServer::~MapMakerServer()
{
}

void MapMakerServer::handle_map_request(const zh::request& request, const el::scope& scope, zh::reply& reply)
{
	string pdbID = request.get_parameter("id", ""s);
	string type = request.get_parameter("type", "density"s);
	string stage = request.get_parameter("stage", "final"s);
	float samplingRate = request.get_parameter("sampling-rate", 1.5f);
	float border = request.get_parameter("border", 5.f);
	bool masked = request.has_parameter("masked");
	
	if (type != "density" and type != "difference" and type != "anomalous")
		throw runtime_error("Invalid type specified, allowed values are density, difference and anomalous");
	
	if (stage != "final" and stage != "besttls" and stage != "0cyc")
		runtime_error("Invalid stage specified, allowed values are final, besttls and 0cyc");
	
	ba::to_lower(pdbID);
	
	if (pdbID.length() != 4)	// sja... TODO!!! fix
		throw runtime_error("not a valid pdbid");
	
	fs::path xyzin = gPDBREDODIR / pdbID.substr(1, 2) / pdbID / (pdbID + "_" + stage + ".cif");
	if (not fs::exists(xyzin))
		xyzin = gPDBREDODIR / pdbID.substr(1, 2) / pdbID / (pdbID + "_" + stage + ".cif.gz");
	if (not fs::exists(xyzin))
		xyzin = gPDBREDODIR / pdbID.substr(1, 2) / pdbID / (pdbID + "_" + stage + ".pdb");
	if (not fs::exists(xyzin))
		xyzin = gPDBREDODIR / pdbID.substr(1, 2) / pdbID / (pdbID + "_" + stage + ".pdb.gz");

	if (not fs::exists(xyzin))
		throw runtime_error("pdb file not found: " + xyzin.string());

	fs::path mtzin = gPDBREDODIR / pdbID.substr(1, 2) / pdbID / (pdbID + "_" + stage + ".mtz");
	if (not fs::exists(mtzin))
		mtzin = gPDBREDODIR / pdbID.substr(1, 2) / pdbID / (pdbID + "_" + stage + ".mtz.gz");
	
	if (not fs::exists(mtzin))
		throw runtime_error("mtz file not found: " + mtzin.string());

	createMap(xyzin, mtzin, type, masked, border, samplingRate, reply);

	string fileName = pdbID + '-' + type;
	if (masked)
		fileName += "-masked";
	fileName += ".map";

    reply.set_header("Content-disposition", "attachment; filename=" + fileName);
}

void MapMakerServer::load_template(const std::string& file, zeep::xml::document& doc)
{
	mrsrc::rsrc rsrc(file);
	if (not rsrc)
		throw runtime_error("missing template");
	
	string data(rsrc.data(), rsrc.size());
	doc.read(data);
}

void MapMakerServer::handle_file(const zh::request& request, const el::scope& scope, zh::reply& reply)
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
	po::options_description visible_options("map-maker " + VERSION_STRING + " options");
	visible_options.add_options()
		("help,h",										"Display help message")
		("hklin",			po::value<string>(),		"Input file (either mtz or cif reflections file)")
		("xyzin",			po::value<string>(),		"Input coordinates file")
		("mapout",			po::value<string>(),		"Output map file")
		("type",			po::value<string>(),		"Map type to generate can be one of density, difference or anomalous")
		
		("sampling-rate",	po::value<float>(),			"Sampling rate, default is 0.75")

		("masked",										"Create a map around the coordinates of the atoms in de model")
		("border",			po::value<float>(),			"Border around the masked region, default is 5")
		
		("no-bulk",										"No bulk ")
		("electron-scattering",							"Use electron scattering factors")
		("aniso-obs",									"Anisotropic scaling of observed")
		("aniso-cal",									"Anisotropic scaling of calculated")

		("server",										"Start as web service")
		("address",				po::value<string>(),	"External address, default is 0.0.0.0")
		("port",				po::value<uint16_t>(),	"Port to listen to, default is 10337")
		("no-daemon,F",									"Do not fork into background")
		("user,u",				po::value<string>(),	"User to run the daemon")
		("logfile",				po::value<string>(),	"Logfile to write to, default /var/log/rama-angles.log")

		("pdb-dir",				po::value<string>(),	"PDB mmCIF data directory")
		("pdb-redo-dir",		po::value<string>(),	"PDB-REDO data directory")
		("ccp4-dir",			po::value<string>(),	"Directory where ccp4 is installed")

		("version",										"Print version")
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		cerr << visible_options << endl;
		exit(0);
	}

	if (vm.count("server") == 0 and
		(vm.count("help") or vm.count("hklin") == 0 or vm.count("xyzin") == 0 or vm.count("mapout") == 0 or vm.count("type") == 0))
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

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
		
		uint16_t port = 10337;
		if (vm.count("port"))
			port = vm["port"].as<uint16_t>();

		string user = "nobody";
		if (vm.count("user"))
			user = vm["user"].as<string>();
		
		bool fork = vm.count("no-daemon") == 0;

		if (fork)
			Daemonize();
		
		string logFile = "/var/log/map-server.log";
		if (vm.count("logfile"))
			logFile = vm["logfile"].as<string>();

		auto serverFactory = []() -> zeep::http::server*
		{
			return new MapMakerServer();
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

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

	mmcif::File f(vm["xyzin"].as<string>());
	mmcif::Structure structure(f);

	mmcif::MapMaker<float> mm;
	
	fs::path hklin = vm["hklin"].as<string>();
	
	if (mmcif::IsMTZFile(hklin))
		mm.loadMTZ(hklin, samplingRate);
	else
	{
		auto aniso = c::MapMaker<float>::as_None;
		if (vm.count("aniso-scaling"))
		{
			if (vm["aniso-scaling"].as<string>() == "observed")
				aniso = c::MapMaker<float>::as_Observed;
			else if (vm["aniso-scaling"].as<string>() == "calculated")
				aniso = c::MapMaker<float>::as_Calculated;
		}

		bool electronScattering = vm.count("electron-scattering") > 0;
		if (not electronScattering)
		{
			auto& exptl = f.data()["exptl"];
			electronScattering = not exptl.empty() and exptl.front()["method"] == "ELECTRON CRYSTALLOGRAPHY";
		}
		
		mm.calculate(hklin, structure, vm.count("no-bulk"), aniso, samplingRate, electronScattering);
	}
	
	string type = vm["type"].as<string>();
	fs::path mapout = vm["mapout"].as<string>();
	fs::ofstream mapFile(mapout, ios_base::binary);
	if (not mapFile.is_open())
		throw runtime_error("Could not open map file " + mapout.string());

	auto range = mm.fb().get().grid_asu();
	if (vm.count("masked"))
	{
		float border = 5;
		if (vm.count("border"))
			border = vm["border"].as<float>();
		range = calculateRange(structure, mm.fb(), border);
	}

	if (type == "density")
		mm.fb().write_masked(mapout, range);
	else if (type == "difference")
		mm.fd().write_masked(mapout, range);
	else if (type == "anomalous")
	{
		if (mm.fa().get().is_null())
			throw runtime_error("Anomalous map not present");
		mm.fa().write_masked(mapout, range);
	}
	else
		throw runtime_error("Unsupported map type " + type);

	return 0;
}
