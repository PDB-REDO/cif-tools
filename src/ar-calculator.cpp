/* 
   Created by: Maarten L. Hekkelman
   Date: 9 oktober 2017
*/

#include "pdb-redo.h"

#include <iostream>
#include <filesystem>s

#include <boost/program_options.hpp>

#include "mmcif/structure.h"

#include "progress.h"

using namespace std;

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

int cif::VERBOSE = 0;

// avoid having to include all the clipper libraries... sigh

namespace clipper
{

int Message::message_level_;
int Message::fatal_level_;
std::ostream* Message::stream_;
	
}

// --------------------------------------------------------------------

template<typename M>
clipper::Mat33<float> normalize(M m)
{
	clipper::Mat33<float> r;
	
	for (int i = 0; i < 3; ++i)
	{
		auto v = clipper::Vec3<float>(m(i, 0), m(i, 1), m(i, 2)).unit();
		r(i, 0) = v[0];
		r(i, 1) = v[1];
		r(i, 2) = v[2];
	}
	
	return r;
}

tuple<float,float> ar_calculate(const mmcif::atom_view& atoms, bool includeWeight)
{
	float If[6] = {};

	for (auto atom: atoms)
	{
		float x, y, z;
		tie(x, y, z) = atom.location();
		
		float w = 1;
		if (includeWeight)
			w = mmcif::atom_type_traits(atom.type()).weight();
		
		If[0] += w * (y * y + z * z);	// 11
		If[1] += w * (x * x + z * z);	// 22
		If[2] += w * (x * x + y * y);	// 33
		If[3] -= w * x * y;				// 12
		If[4] -= w * x * z;				// 13
		If[5] -= w * y * z;				// 23
	}

	auto m = clipper::Mat33<float>(
		clipper::Mat33sym<float>(If[0], If[1], If[2], If[3], If[4], If[5])
	);

	clipper::Matrix<float> M(3, 3);

	auto mn = normalize(m);
	
	M(0, 0) = mn(0, 0);
	M(1, 1) = mn(1, 1);
	M(2, 2) = mn(2, 2);
	M(0, 1) = M(1, 0) = mn(0, 1);
	M(0, 2) = M(2, 0) = mn(0, 2);
	M(1, 2) = M(2, 1) = mn(1, 2);
	
	vector<float> eigen_values = M.eigen();
				
	sort(eigen_values.begin(), eigen_values.end());

	float axiality = 2 * eigen_values[2] - (eigen_values[0] + eigen_values[1]);
//	float rhombicity = eigen_values[2] + eigen_values[1];
	float rhombicity = eigen_values[1] - eigen_values[0];
	
	return make_tuple(axiality, rhombicity);
}

// --------------------------------------------------------------------

void ar_calculator(fs::path file)
{
	if (cif::VERBOSE)
		cout << string(get_terminal_width(), '=') << endl
			 << "Processing file: " << file << endl;
	
	c::file mmcif(file);
	c::structure prot(mmcif, 1);

	float axiality, axiality_weighted, rhombicity, rhombicity_weighted;
	tie(axiality, rhombicity) = ar_calculate(prot.atoms(), false);
	tie(axiality_weighted, rhombicity_weighted) = ar_calculate(prot.atoms(), true);
	
	cout << file << '\t' << axiality << '\t' << rhombicity << '\t' << axiality_weighted << '\t' << rhombicity_weighted << endl;
}

// --------------------------------------------------------------------

int main(int argc, char* argv[])
{
	po::options_description desc("ar_calculator " + VERSION_STRING + " options");
	desc.add_options()
		("input,i",		po::value<string>(),	"Input file or directory")
		("recursive,r",							"Process recusively on directory")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::positional_options_description p;
	p.add("input", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << desc << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	fs::path file = vm["input"].as<string>();

	cout << "filename" << '\t' << "axiality" << '\t' << "rhombicity" << '\t' << "axiality(w)" << '\t' << "rhombicity(w)" << endl
		 << string(get_terminal_width(), '-') << endl;
	
	if (fs::is_regular(file))
		ar_calculator(file);
	else if (fs::is_directory(file))
	{
		if (vm.count("recursive"))
		{
			for (auto i = fs::recursive_directory_iterator(file); i != fs::recursive_directory_iterator(); ++i)
			{
				if (fs::is_regular(*i))
					ar_calculator(*i);
			}
		}
		else
		{
			for (auto i = fs::directory_iterator(file); i != fs::directory_iterator(); ++i)
			{
				if (fs::is_regular(*i))
					ar_calculator(*i);
			}
		}
	}
	else
	{
		cerr << "Could not open " << file << endl;
		exit(1);
	}

	return 0;
}
