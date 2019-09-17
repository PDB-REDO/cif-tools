/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 03 oktober, 2016
*/

// NOTE: TODO the OXT atoms are ignored at multiple locations


#include "pdb-redo.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <numeric>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#include <zeep/xml/document.hpp>
#include <zeep/xml/serialize.hpp>
#include <zeep/xml/writer.hpp>

#include "cif++/Cif++.h"
#include "cif++/Structure.h"
#include "cif++/Compound.h"
//#include "cif++/Edia.h"
#include "cif++/CifUtils.h"
#include "cif++/DistanceMap.h"
#include "cif++/BondMap.h"
#include "cif++/MapMaker.h"

#include "HBondTraits.h"

#include "cif++/mrsrc.h"
#include "svm++.h"

using namespace std;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace zx = zeep::xml;
namespace c = mmcif;

// --------------------------------------------------------------------
// global constants, or magic numbers (may need to change these to variables?)

const static float
	kMinSigmaToNormalise = 0.01f,
	kRaiseMapSigmas = 0.f,
	kAtomBumpDensity = 2.0f,
	kAtomBumpThreshold = 0.05f,
	kAtomRadiusIncreaseFactor = 1.2f,
	kBFactorToAtomRadiusExtensionFactor = 79.f,
	kInflationFactor = 1;

// --------------------------------------------------------------------

enum WaterFeaturesEnum {
	kFeatureDensity,
	kFeatureDensityIntegrated,
	kFeatureInterpolatedWaterRadius,
	kFeatureEDIA,
	kFeatureDensitySDAtRadius,
	kFeatureAxiality,
	kFeatureRhombicity,
//	kFeatureNrOfBumps,
	kFeatureNrOfHBonds,
//	kFeatureSumOfInvHBondLengths,
	kFeatureEstimatedHBondEnergy,
	
	kNrOfFeatures
};

typedef array<float,kNrOfFeatures> FeatureScoreArray;

// --------------------------------------------------------------------

struct CentrifugeParameters
{
	float					uniformAtomRadius;
	float					uniformBFactor;
	float					minWaterDistance;
	float					maxWaterDistance;
	float					waterSphereRadius;
	float					waterDensityRadiusCutOff;
	vector<string>			mtzColumns;
	float					featureScales[kNrOfFeatures];
	float					featureOffsets[kNrOfFeatures];
	
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		ar & ZEEP_ELEMENT_NAME_VALUE(uniformAtomRadius)
		   & ZEEP_ELEMENT_NAME_VALUE(uniformBFactor)
		   & ZEEP_ELEMENT_NAME_VALUE(minWaterDistance)
		   & ZEEP_ELEMENT_NAME_VALUE(maxWaterDistance)
		   & ZEEP_ELEMENT_NAME_VALUE(waterSphereRadius)
		   & ZEEP_ELEMENT_NAME_VALUE(waterDensityRadiusCutOff)
		   & ZEEP_ELEMENT_NAME_VALUE(mtzColumns);
	}
};

const static CentrifugeParameters
	kDefaultParameters = {
		0,					// uniformAtomRadius;
		20,					// uniformBFactor;
		2.2f,				// minWaterDistance;
		3.2f,				// maxWaterDistance;
		0,					// waterSphereRadius
		0.15f,				// waterDensityRadiusCutOff;
		{ "FWT", "PHWT" }	// mtzColumns;
	};

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

// --------------------------------------------------------------------

void CalculateIntegrationBox(const clipper::Xmap<float>& xmap, const vector<c::Atom>& atoms,
	clipper::Coord_grid& gridMin, clipper::Coord_grid& gridMax,
	float inflationFactor, float uniformRadius, float uniformBFactor)
{
	// calculate min and max orthogonal coordinates first, based on atom position, radius and bfactor
	
	clipper::Coord_orth oMin, oMax;
	
	oMin[0] = oMin[1] = oMin[2] = numeric_limits<float>::max();
	oMax[0] = oMax[1] = oMax[2] = numeric_limits<float>::lowest();
	
	for (auto& atom: atoms)
	{
		float x, y, z;
		tie(x, y, z) = atom.location();
		
		if (isnan(x) or isnan(y) or isnan(z))
			continue;

		float radius = uniformRadius;
		if (radius == 0)
			radius = c::AtomTypeTraits(atom.type()).radius();
		radius = sqrt(radius * radius + uniformBFactor / kBFactorToAtomRadiusExtensionFactor) * inflationFactor;
		
		if (oMin[0] > x - radius * 2)
			oMin[0] = x - radius * 2;
		if (oMax[0] < x + radius * 2)
			oMax[0] = x + radius * 2;

		if (oMin[1] > y - radius * 2)
			oMin[1] = y - radius * 2;
		if (oMax[1] < y + radius * 2)
			oMax[1] = y + radius * 2;

		if (oMin[2] > z - radius * 2)
			oMin[2] = z - radius * 2;
		if (oMax[2] < z + radius * 2)
			oMax[2] = z + radius * 2;
	}
	
	clipper::Coord_frac fMin = oMin.coord_frac(xmap.cell());
	clipper::Coord_frac fMax = oMax.coord_frac(xmap.cell());

	clipper::Coord_map mMin = fMin.coord_map(xmap.grid_sampling());
	clipper::Coord_map mMax = fMax.coord_map(xmap.grid_sampling());

	gridMin = mMin.floor();
	gridMax = mMax.ceil();
}

void GetGridOrthoCoordinates(clipper::Xmap<float>& xmap,
	clipper::Coord_grid& gridMin, clipper::Coord_grid& gridMax, 
	vector<clipper::Coord_orth>& coords)
{
	auto dim = gridMax - gridMin;
	size_t numOfGridPoints = (dim[0] + 1) * (dim[1] + 1) * (dim[2] + 1);

	coords.clear();
	coords.reserve(numOfGridPoints);
	
// TODO: vervangen door iterateGrid
#pragma warning("vervangen door iterateGrid")
	
	auto i0 = clipper::Xmap_base::Map_reference_coord(xmap, gridMin);
	for (auto iu = i0; iu.coord().u() <= gridMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gridMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gridMax[2]; iw.next_w())
				coords.push_back(iw.coord_orth());
}

//void RemoveAtoms(clipper::Xmap<float>& xmap, const vector<c::Atom>& atoms,
//	float inflationFactor, float uniformRadius = 0, float uniformBFactor = 0)
//{
//	clipper::Coord_grid gridMin, gridMax;
//	
//	CalculateIntegrationBox(xmap, atoms, gridMin, gridMax, inflationFactor, uniformRadius, uniformBFactor);
//
//	vector<clipper::Coord_orth> coordinates;
//	GetGridOrthoCoordinates(xmap, gridMin, gridMax, coordinates);
//
//	vector<float> modDensity(coordinates.size(), 0.0f);
//	
//	for (const c::Atom& atom: atoms)
//	{
//		float radius = uniformRadius;
//		if (radius == 0)
//			radius = c::AtomTypeTraits(atom.type()).radius();
//		radius = sqrt(radius * radius + uniformBFactor / kBFactorToAtomRadiusExtensionFactor) * inflationFactor;
//
//		size_t i = 0;
//		for (auto c: coordinates)
//		{
//			float x, y, z;
//			tie(x, y, z) = atom.location();
//
//			float d2 =
//				(x - c[0]) * (x - c[0]) +
//				(y - c[1]) * (y - c[1]) +
//				(z - c[2]) * (z - c[2]);
//			
//			d2 /= radius * radius;
//			
//			float d8 = 0.3025f * d2;
//			d8 = d8 * d8 * d8 * d8;
//			
//			modDensity[i] += exp(-d2 - d8);
//			++i;
//		}
//	}
//	
//	auto i0 = clipper::Xmap_base::Map_reference_coord(xmap, gridMin);
//	size_t i = 0;
//	
//	for (auto iu = i0; iu.coord().u() <= gridMax[0]; iu.next_u())
//	{
//		for (auto iv = iu; iv.coord().v() <= gridMax[1]; iv.next_v())
//		{
//			for (auto iw = iv; iw.coord().w() <= gridMax[2]; iw.next_w())
//			{
//				assert(i < modDensity.size());
//				
//				modDensity[i] *= -kAtomBumpDensity;
//
//				if (modDensity[i] < -kAtomBumpThreshold and modDensity[i] < xmap[iw])
//					xmap[iw] = modDensity[i];
//				++i;
//			}
//		}
//	}
//}

// --------------------------------------------------------------------

float CalculateDensity(const clipper::Xmap<float>& xmap, const vector<c::Atom>& atoms,
	float inflationFactor, float uniformRadius, float uniformBFactor)
{
	clipper::Coord_grid gridMin, gridMax;
	
	CalculateIntegrationBox(xmap, atoms, gridMin, gridMax, inflationFactor, uniformRadius, uniformBFactor);
	
	double sumDensity = 0, sumCorrelation = 0;
	
	for (const c::Atom& atom: atoms)
	{
		float radius = uniformRadius;
		if (radius == 0)
			radius = c::AtomTypeTraits(atom.type()).radius();
		
		radius = sqrt(radius * radius + uniformBFactor / kBFactorToAtomRadiusExtensionFactor);

// TODO: vervangen door iterateGrid
#pragma warning("vervangen door iterateGrid")
	
		auto i0 = clipper::Xmap_base::Map_reference_coord(xmap, gridMin);
		for (auto iu = i0; iu.coord().u() <= gridMax[0]; iu.next_u())
			for (auto iv = iu; iv.coord().v() <= gridMax[1]; iv.next_v())
				for (auto iw = iv; iw.coord().w() <= gridMax[2]; iw.next_w())
				{
					auto c = iw.coord_orth();
					
					float x, y, z;
					tie(x, y, z) = atom.location();

					float d2 =
						(x - c[0]) * (x - c[0]) +
						(y - c[1]) * (y - c[1]) +
						(z - c[2]) * (z - c[2]);
					
					d2 /= radius * radius;
					
					float d8 = 0.3025f * d2;
					d8 = d8 * d8 * d8 * d8;
					
					float v = exp(-d2 - d8);
					
					sumDensity += v;
					sumCorrelation += v * xmap[iw];
				}
	}
	
	float result;
	
	if (sumDensity > 0)
		result = static_cast<float>(sumCorrelation / sumDensity);
	else
		result = 0;

	return result;
}

// --------------------------------------------------------------------

float InterpolateDensityRadius(const clipper::Xmap<float>& xmap, const c::Atom& atom, float density)
{
	static const c::SphericalDots_50 dots;

	struct {
		float	radius;
		double	avg;
		double	sd;
	} samples[] = {
		 { 0.15f },
		 { 0.3f },
		 { 0.45f },
		 { 0.6f },
		 { 0.75f },
		 { 0.9f },
		 { 1.05f },
		 { 1.2f },
		 { 1.35f },
		 { 1.5f },
		 { 1.65f }
	};
	
	for (auto& sample: samples)
	{
		float radius = sample.radius;

		vector<double> d;
		d.reserve(dots.size());
		double sum = 0;
		
		for (auto dot: dots)
		{
			dot *= radius;
			
			c::Point p = atom.location() + dot;
	
			clipper::Coord_orth co = p;
			clipper::Coord_frac cf = co.coord_frac(xmap.cell());
			
			float dp = xmap.interp<clipper::Interp_cubic>(cf);
			
			d.push_back(dp);
			sum += dp;
		}
		
		sample.avg = sum / dots.size();
		double sum2 = accumulate(d.begin(), d.end(),
			0., [=](double s, double d) -> double
			{
				return s + pow(d - sample.avg, 2);
			}
		);
		
		sample.sd = sqrt(sum2 / dots.size());
		
		if (VERBOSE > 1)
			cout << "For r = " << radius << " avg: " << sample.avg << " and sd " << sample.sd << endl;
	}
	
	if (VERBOSE > 1)
		cout << endl;

	// Use linear interpolation for now

	clipper::Coord_orth co = atom.location();
	clipper::Coord_frac cf = co.coord_frac(xmap.cell());
	
	float result = 0;
	
	float r = 0;
	float d = xmap.interp<clipper::Interp_cubic>(cf);
	
	for (auto sample: samples)
	{
		if (density <= d and density >= sample.avg)
		{
			float f = (density - d) / (sample.avg - d);
			result = r + f * (sample.radius - r);
			break;
		}
		
		r = sample.radius;
		d = sample.avg;
	}
	
	return result;
}

// --------------------------------------------------------------------

float DensityStandardDeviationAtRadius(const clipper::Xmap<float>& xmap, const c::Atom& atom, float radius)
{
	static const c::SphericalDots_50 dots;

	vector<double> d;
	d.reserve(dots.size());
	double sum = 0;
		
	for (auto dot: dots)
	{
		dot *= radius;
		
		c::Point p = atom.location() + dot;

		clipper::Coord_orth co = p;
		clipper::Coord_frac cf = co.coord_frac(xmap.cell());
		
		float dp = xmap.interp<clipper::Interp_cubic>(cf);
		
		d.push_back(dp);
		sum += dp;
	}

	auto avg = sum / dots.size();
	double sum2 = accumulate(d.begin(), d.end(),
		0., [=](double s, double d) -> double
		{
			return s + pow(d - avg, 2);
		}
	);
		
	auto sd = sqrt(sum2 / dots.size());
	
	return static_cast<float>(sd);
}

//// --------------------------------------------------------------------
//
//clipper::Mat33<float> CalculateInertiaTensor(const clipper::Xmap<float>& xmap, const c::Atom& atom, float radius)
//{
//	static const c::SphericalDots_50 dots;
//	float If[6] = {};
//	
//	for (auto dot: dots)
//	{
//		dot *= radius;
//		
//		c::Point p = atom.location() + dot;
//
//		clipper::Coord_orth co = p;
//		clipper::Coord_frac cf = co.coord_frac(xmap.cell());
//		float dp = xmap.interp<clipper::Interp_cubic>(cf);
//		
//		float x, y, z;
//		tie(x, y, z) = p;
//		
//		If[0] += dp * (y * y + z * z);	// 11
//		If[1] += dp * (x * x + z * z);	// 22
//		If[2] += dp * (x * x + y * y);	// 33
//		If[3] -= dp * x * y;			// 12
//		If[4] -= dp * x * z;			// 13
//		If[5] -= dp * y * z;			// 23
//	}
//
//	auto result = clipper::Mat33<float>(
//		clipper::Mat33sym<float>(If[0], If[1], If[2], If[3], If[4], If[5])
//	);
//
//	return result;
//}
//

// --------------------------------------------------------------------

clipper::Mat33<float> CalculateInertiaTensor(const clipper::Xmap<float>& xmap,
	const c::Atom& atom, float radius, float densityCutoff)
{
	float If[6] = {};
	auto c= atom.location();
	
	float x, y, z;
	tie(x, y, z) = c;

	// calculate min and max orthogonal coordinates first, based on atom position, radius and bfactor
	clipper::Coord_orth oMin = { x - radius * 2, y - radius * 2, z - radius * 2},
						oMax = { x + radius * 2, y + radius * 2, z + radius * 2};
	
	clipper::Coord_frac fMin = oMin.coord_frac(xmap.cell());
	clipper::Coord_frac fMax = oMax.coord_frac(xmap.cell());

	clipper::Coord_map mMin = fMin.coord_map(xmap.grid_sampling());
	clipper::Coord_map mMax = fMax.coord_map(xmap.grid_sampling());

	clipper::Coord_grid gridMin = mMin.floor();
	clipper::Coord_grid gridMax = mMax.ceil();

// TODO: vervangen door iterateGrid
#pragma warning("vervangen door iterateGrid")

	auto i0 = clipper::Xmap_base::Map_reference_coord(xmap, gridMin);
	for (auto iu = i0; iu.coord().u() <= gridMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gridMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gridMax[2]; iw.next_w())
			{
				float dp = xmap[iw];
				
				if (dp < densityCutoff)	
					continue;
					
				tie(x, y, z) = c::Point(iw.coord_orth());
				
				If[0] += dp * (y * y + z * z);	// 11
				If[1] += dp * (x * x + z * z);	// 22
				If[2] += dp * (x * x + y * y);	// 33
				If[3] -= dp * x * y;			// 12
				If[4] -= dp * x * z;			// 13
				If[5] -= dp * y * z;			// 23
			}
		

	auto result = clipper::Mat33<float>(
		clipper::Mat33sym<float>(If[0], If[1], If[2], If[3], If[4], If[5])
	);

	return result;
}

// --------------------------------------------------------------------
// same as above, however this time we use all grid points within the radius

clipper::Mat33<float> CalculateInertiaTensor2(const clipper::Xmap<float>& xmap, const c::Atom& atom,
	float waterRadius, float inflationFactor, float uniformRadius, float uniformBFactor)
{
	clipper::Coord_grid gridMin, gridMax;
	
	CalculateIntegrationBox(xmap, c::AtomView({atom}), gridMin, gridMax, inflationFactor, uniformRadius, uniformBFactor);
	
//	double sumDensity = 0, sumCorrelation = 0;
	
	float If[6] = {};

// TODO: vervangen door iterateGrid
#pragma warning("vervangen door iterateGrid")

	auto i0 = clipper::Xmap_base::Map_reference_coord(xmap, gridMin);
	for (auto iu = i0; iu.coord().u() <= gridMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gridMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gridMax[2]; iw.next_w())
			{
				auto c = iw.coord_orth();
				
				float distance = c::Distance(atom.location(), c::Point{c});
				
				if (distance > waterRadius)
					continue;
				
				float dp = xmap[iw];
				
				If[0] += dp * (c[1] * c[1] + c[2] * c[2]);	// 11
				If[1] += dp * (c[0] * c[0] + c[2] * c[2]);	// 22
				If[2] += dp * (c[0] * c[0] + c[1] * c[1]);	// 33
				If[3] -= dp * c[0] * c[1];					// 12
				If[4] -= dp * c[0] * c[2];					// 13
				If[5] -= dp * c[1] * c[2];					// 23
			}

	auto result = clipper::Mat33<float>(
		clipper::Mat33sym<float>(If[0], If[1], If[2], If[3], If[4], If[5])
	);

	return result;
}

// --------------------------------------------------------------------

tuple<float,float> AxialityAndRhombicity(const clipper::Mat33<float>& m)
{
	clipper::Matrix<float> M(3, 3);

	auto mn = normalize(m);
	
	M(0, 0) = mn(0, 0);
	M(1, 1) = mn(1, 1);
	M(2, 2) = mn(2, 2);
	M(0, 1) = M(1, 0) = mn(0, 1);
	M(0, 2) = M(2, 0) = mn(0, 2);
	M(1, 2) = M(2, 1) = mn(1, 2);
	
	vector<float> eigenValues = M.eigen();
				
	sort(eigenValues.begin(), eigenValues.end());

	float axiality = 2 * eigenValues[2] - (eigenValues[0] + eigenValues[1]);
	float rhombicity = eigenValues[1] - eigenValues[0];
	
	return make_tuple(axiality, rhombicity);
}

// --------------------------------------------------------------------

void SaveMap(clipper::Xmap<float>& xmap, fs::path f)
{
	clipper::CCP4MAPfile tmpMapOutFile;
	tmpMapOutFile.open_write(f.string());
	tmpMapOutFile.export_xmap(xmap);
	tmpMapOutFile.close_write();
}

// --------------------------------------------------------------------

struct score
{
	string						name;
	c::Point					loc;
	float						density, densityIntegrated;
	float						interpolatedWaterRadius;
	clipper::Mat33<float>		inertiaTensor;
	float						axiality, rhombicity;
	vector<pair<string,float>>	bumps, hbondpairs;
	float						edia;
	
	static string columnHeaders()
	{
		stringstream s;
		s << "Name" << '\t'
		  << "Loc" << '\t'
		  << "Density" << '\t'
		  << "Density integrated" << '\t'
		  << "EDIA" << '\t'
		  << "Interpolated water radius" << '\t'
		  << "Axiality" << '\t'
		  << "Rhombicity" << '\t'
		  << "Bumps" << '\t'
		  << "HBonds";
		return s.str();
	}
};

ostream& operator<<(ostream& os, const score& sc)
{
	os << sc.name << '\t'
	   << sc.loc << '\t'
	   << sc.density << '\t'
	   << sc.densityIntegrated << '\t'
	   << sc.edia << '\t';
	
	if (sc.interpolatedWaterRadius > 0)
		os << sc.interpolatedWaterRadius << '\t'
		   << sc.axiality << '\t'
		   << sc.rhombicity << '\t';
	else
		os << '\t'
		   << '\t'
		   << '\t';

	os << sc.bumps.size() << '\t'
	   << sc.hbondpairs.size();		
	
	return os;
}

// --------------------------------------------------------------------

FeatureScoreArray CalculateScoresForWater(const c::Atom& water, const c::Structure& prot,
	const clipper::Xmap<float>& xmap, const c::DistanceMap& distances, const c::BondMap& bonds,
	float resolution, float meanDensity, float rmsDensity,
	NotAHBondSet& notAHBond, const CentrifugeParameters& params)
{
	if (VERBOSE)
		cout << string(cif::get_terminal_width(), '-') << endl;

//	if (not water.isWater())
//		throw runtime_error("Invalid water (" + water.mChainID + to_string(water.mResSeq) + "), not an oxygen or resname is not HOH");
	
	clipper::Coord_orth co = water.location();
	clipper::Coord_frac cf = co.coord_frac(xmap.cell());

	FeatureScoreArray sc;
	sc.fill(0);
	
	sc[kFeatureDensity] = xmap.interp<clipper::Interp_cubic>(cf);
	sc[kFeatureDensityIntegrated] = CalculateDensity(xmap, vector<c::Atom>({ water }),
		kInflationFactor, params.uniformAtomRadius, params.uniformBFactor);
	
//	sc[kFeatureEDIA] = CalculateEDIA(water, xmap, resolution, meanDensity, rmsDensity, distances, bonds);
	
	sc[kFeatureInterpolatedWaterRadius] = InterpolateDensityRadius(xmap, water, params.waterDensityRadiusCutOff);

	float radius = params.waterSphereRadius;
	if (radius == 0)
		radius = sc[kFeatureInterpolatedWaterRadius];

	if (radius > 0)
	{
		if (sc[kFeatureInterpolatedWaterRadius] == 0)
			sc[kFeatureInterpolatedWaterRadius] = nan("");
		
		sc[kFeatureDensitySDAtRadius] = DensityStandardDeviationAtRadius(xmap, water, radius);
	
		auto inertiaTensor = CalculateInertiaTensor(xmap, water, radius, params.waterDensityRadiusCutOff);
		tie(sc[kFeatureAxiality], sc[kFeatureRhombicity]) = AxialityAndRhombicity(inertiaTensor);
	}
	else
	{
		sc[kFeatureInterpolatedWaterRadius] = nan("");
		sc[kFeatureDensitySDAtRadius] = nan("");
		sc[kFeatureAxiality] = nan("");
		sc[kFeatureRhombicity] = nan("");
	}

	for (auto& atom: prot.atoms())
	{
		// skip the current water atom
		if (atom.location() == water.location())
			continue;

		float distance = distances(atom, water);

		if (distance < params.minWaterDistance)
		{
			if (VERBOSE)
			{
				cout << "Bumps with atom " << atom.id() << endl;
				if (VERBOSE > 1)
					cout << "atom: " << atom.location() << endl
						 << "water: " << water.location() << endl
						 << "distance = " << distance << endl;
			}
					
//			sc[kFeatureNrOfBumps] += 1;
		}
		else if (distance < params.maxWaterDistance and (atom.type() == c::O or atom.type() == c::N))
		{
			if (notAHBond(atom))
			{
				if (VERBOSE > 1)
					cout << "Atom " << atom.id() << " is not a HBond donor or acceptor" << endl;
			}
			else
			{
				if (VERBOSE > 1)
					cout << "HBond with " << atom.id() << " at " << distance << "" << endl;
	
				sc[kFeatureNrOfHBonds] += 1;
//				sc[kFeatureSumOfInvHBondLengths] += 1 / distance;

				double energyEstimate = 5 * pow(2 / distance, 12) - 6 * pow(2 / distance, 10);
				
				sc[kFeatureEstimatedHBondEnergy] += energyEstimate;
			}
		}
	};

	if (VERBOSE)
	{
		cout << "density at point: " << sc[kFeatureDensity] << endl
			 << "density integrated: " << sc[kFeatureDensityIntegrated] << endl
			 << "EDIA: " << sc[kFeatureEDIA] << endl
			 << "interpolated water radius: " << sc[kFeatureInterpolatedWaterRadius] << endl
			 << "stdev density at radius: " << sc[kFeatureDensitySDAtRadius] << endl
			 << "axiality: " << sc[kFeatureAxiality] << endl
			 << "rhombicity: " << sc[kFeatureRhombicity] << endl
//			 << "Bump: " << sc[kFeatureNrOfBumps] << endl
			 << "HBonds: " << sc[kFeatureNrOfHBonds] << endl
//			 << "Sum of Inv HBonds: " << sc[kFeatureSumOfInvHBondLengths] << endl;
 			 << "HBond energy: " << sc[kFeatureEstimatedHBondEnergy] << endl;
	}

	return sc;
}

// --------------------------------------------------------------------

int centrifugeLearn(int argc, char* argv[])
{
	po::options_description desc("centrifuge-learn " + VERSION + " options");
	desc.add_options()
		("output,o",	po::value<string>(),	"Output file for the calculated model")
		("input,i",		po::value<vector<string>>(),
												"Input files, can be specified multiple times")
												
		("force,f",								"Force overwriting files")
//		("report,r",	po::value<string>(),	"Write the statics report to this file, use 'stdout' to output to screen")

		("dump-features",
						po::value<string>(),	"Dump the final feature table to this file")

		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)")

		// specific options for centrifuge

		("remove-residues",						"Remove residues from map")

		("uniform-atom-radius",
						po::value<float>(),		"Uniform Atom Radius to use for all atoms, default is to use single covalent bond length")

		("uniform-bfactor",
						po::value<float>(),		"Uniform B Factor to use for all atoms, default is 20")

//		("save-fit-map",
//						po::value<string>(),	"Save fit map after residue removal to this file")

		("min-water-distance",
						po::value<float>(),		"Minimum distance between a water atom and any atom of the protein structure, default is 2.0 Å")

		("max-water-distance",
						po::value<float>(),		"Maximum distance between a water atom and any atom of the protein structure, default is 3.5 Å")

		("water-sphere-radius",
						po::value<float>(),		"The radius of the sphere around water to calculate moments of intertia, default is to use calculated value.")

		("water-density-radius-cutoff",
						po::value<float>(),		"The density cutoff to use in calculating moments of intertia, default is 0.75 Å")

		("mtz-columns",
						po::value<string>(),	"Comma separated list of MTZ Column labels to load (default value is \"FWT,PHWT\"")

		("include-disputed",					"Include all waters in model, not only consensus waters")

		("grid-search",
												"Do a grid search for the best C and gamma parameters (expensive!)")

		("match-radius", po::value<float>(),	"Match waters that are at most this much relocated in curated files, default is 0.5 Å")

		("sampling-rate",po::value<float>(),	"Sampling rate")

		;

	po::positional_options_description p;
//	p.add("output", 1);
	p.add("input", -1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("output") == 0 or vm.count("input") == 0)
	{
		cerr << desc << endl;
		exit(1);
	}
	
	if (fs::exists(vm["output"].as<string>()) and vm.count("force") == 0)
	{
		cerr << "cowardly refusing to overwrite " << vm["output"].as<string>() << endl;
		exit(1);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();
	
	// collect parameters
	CentrifugeParameters params = kDefaultParameters;
	if (vm.count("uniform-atom-radius"))
		params.uniformAtomRadius = vm["uniform-atom-radius"].as<float>();
	if (vm.count("uniform-bfactor"))
		params.uniformBFactor = vm["uniform-bfactor"].as<float>();
	if (vm.count("min-water-distance"))
	       params.minWaterDistance = vm["min-water-distance"].as<float>();
	if (vm.count("max-water-distance"))
	       params.maxWaterDistance = vm["max-water-distance"].as<float>();
	if (vm.count("water-sphere-radius"))
		params.waterSphereRadius = vm["water-sphere-radius"].as<float>();
	if (vm.count("water-density-radius-cutoff"))
		params.waterDensityRadiusCutOff = vm["water-density-radius-cutoff"].as<float>();
	if (vm.count("mtz-columns"))
		ba::split(params.mtzColumns, vm["mtz-columns"].as<string>(), ba::is_any_of(","));
		
	bool includeDisputed = vm.count("include-disputed") > 0;

	float closeEnoughToBeSame = 0.5f;
	if (vm.count("match-radius"))
		closeEnoughToBeSame = vm["match-radius"].as<float>();
	float closeEnoughToBeSameSq = pow(closeEnoughToBeSame, 2);

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

//	typedef struct { string chainID; int16_t resSeqNr; int label; } WaterInfo;
//	vector<WaterInfo> waters;

	typedef struct
	{
		fs::path					file;
		string						curator;
		shared_ptr<c::File>			mmcif;
		shared_ptr<c::Structure>	prot;
		size_t						nrOfWaters;
	} StructureInfo;

	map<string,vector<StructureInfo>> structures;

	const regex rx("(.{4})_besttls(?:_dry_(.*))?\\.pdb");
	
	// First collect the PDB structures and their curated files
	// iterate over the input files to collect structures
	for (fs::path input: vm["input"].as<vector<string>>())
	{
		smatch m;
		string filename = input.filename().string();
		
		if (VERBOSE)
			cerr << "reading " << filename << endl;
		
		if (not regex_match(filename, m, rx))
			throw runtime_error("Filenames must match (.{4})_besttls(.*)\\.pdb");

		string pdbID = m[1].str();
		string curator = m[2].str();

		StructureInfo si = { input, curator };
		
		si.mmcif = make_shared<c::File>(input);
		si.prot = make_shared<c::Structure>(*si.mmcif, 1);
		si.nrOfWaters = si.prot->waters().size();
		
		structures[pdbID].push_back(si);
	}
	
	// Create a NotAHBondSet dataset
	auto notAHBond = NotAHBondSet::Create();
	
	// process the files, collect data
	svm::Matrix data;
	vector<svm::label_type> labels;

	// Statistics
	typedef struct {
		string							structure;
		uint32_t							waters;
		uint32_t							rejected;
		uint32_t							accepted;
		uint32_t							disputed;
		map<string,uint32_t>				acceptedBy;
	} StructureStatistics;
	vector<StructureStatistics> stats;

	// Dump features, write out a dump of all features to allow analysis
	
	ofstream dump;
	if (vm.count("dump-features"))
	{
		dump.open(vm["dump-features"].as<string>());
		if (not dump.is_open())
			throw runtime_error("Could not write to dump-features file");
	}

	for (auto& s: structures)
	{
		cout << string(cif::get_terminal_width(), '=') << endl;

		auto& proteins = s.second;
		
		if (proteins.size() == 1)
		{
			cout << "Skipping " << s.first << " since there is no curated file yet" << endl;
			continue;
		}
		
		// sort by curator name, since there is no curator for the 'main' it will end up at the
		// head of the list
		sort(proteins.begin(), proteins.end(), [](const StructureInfo& a, const StructureInfo& b) -> bool
			{ return a.curator < b.curator; });
		
		cout << "Chose " << proteins.front().file.filename() << " as main structure for " << s.first << endl
			 << "Nr of curated files: " << (proteins.size() - 1) << endl;
		
		// Note, id's of atom need not be the same across the various files.
		// The approach is therefore to lookup waters by residue sequence number.
		// If that fails, we fall back to search by location.
		map<string,uint32_t> curatedCount;
		struct WaterLoc {
			string		id;
			string		asymId;
			int			resNum;
			c::Point	loc;
		};
		vector<WaterLoc> waterLocations;

		for (auto& w: proteins.front().prot->waters())
		{
			curatedCount[w.id()] = 1;
			waterLocations.push_back(
			{
				w.id(),
				w.property<string>("auth_asym_id"),
				w.property<int>("auth_seq_id"),
				w.location()
			});
		}

		StructureStatistics stat = { s.first };
		
		for (size_t cpi = 1; cpi < proteins.size(); ++cpi)
		{
			auto& cp = proteins[cpi];
			
			uint32_t n = 0;
			for (auto& w: cp.prot->waters())
			{
				string asymId = w.property<string>("auth_asym_id");
				int seqId = w.property<int>("auth_seq_id");
				
				auto i = find_if(waterLocations.begin(), waterLocations.end(), [=](auto& wl) -> bool
					{
						bool result = false;
						if (wl.asymId == asymId and wl.resNum == seqId)
						{
							float d = DistanceSquared(w.location(), wl.loc);
							result = d <= closeEnoughToBeSameSq;
						}
						return result;
					});
				
//				if (i == waterLocations.end())
//				{
//					cerr << "Locating water by id failed for " << asymId << ':' << seqId << endl;
//					
//					i = find_if(waterLocations.begin(), waterLocations.end(), [&](auto& wl) -> bool
//						{
//							float d = DistanceSquared(w.location(), wl.loc);
//							return d <= closeEnoughToBeSameSq;
//						});
//				}

				if (i == waterLocations.end())
					cerr << "unknown (new?) water " << asymId << ':' << seqId << endl;
				else
				{
					if (VERBOSE)
					{
						float d = DistanceSquared(w.location(), i->loc);
						if (d > 0)
							cerr << "Water " << asymId << ':' << seqId << " drifted by " << sqrt(d) << " in " << cp.file << endl;
					}

					curatedCount[i->id] += 1;
					++n;
				}
			}
			
			stat.acceptedBy[cp.curator] = n;
		}
		
		for (auto curated: curatedCount)
		{
			if (curated.second == 1)
				stat.rejected += 1;
			if (curated.second == proteins.size())
				stat.accepted += 1;
		}
		
		stat.disputed = curatedCount.size() - stat.rejected - stat.accepted;
		stats.push_back(stat);
		
		cout << "Nr of waters: " << curatedCount.size() << endl
			 << "Nr of real waters: " << stat.accepted << endl
			 << "Nr of not waters: " << stat.rejected << endl
			 << "Nr of disputed waters: " << stat.disputed << endl;
		
		// Do the actual calculations
		c::Structure& prot = *proteins.front().prot.get();

		// Read the xmap from the corresponding mtz file
		fs::path file = proteins.front().file;
		fs::path mtzFile = file.parent_path() / (file.stem().string() + ".mtz");

		c::MapMaker<float> mm;
		mm.loadMTZ(mtzFile, samplingRate);
	
		const clipper::Xmap<float>& xmap = mm.fb();
		
		float resolution = mm.resHigh(), meanDensity = mm.fb().meanDensity(), rmsDensity = mm.fb().rmsDensity();

		// Calculate the distance map
		auto& spacegroup = xmap.spacegroup();
		auto& cell = xmap.cell();
		
		c::DistanceMap distances(prot, spacegroup, cell, 3.5f);
		c::BondMap bonds(prot);
		
		cif::Progress progress(curatedCount.size(), "Processing waters for " + s.first);
		
		for (auto w: curatedCount)
		{
			auto id = w.first;
			uint32_t cnt = w.second;
			
			
			progress.message(id);
			
			c::Atom water = prot.getAtomById(id);
			
			if (not water.comp().isWater())
				throw runtime_error("This is not a water in " + file.string() + " id: " + id);

			if (VERBOSE)
				cout << string(cif::get_terminal_width(), '=') << endl
					 << "Learning from " << id << " @ " << water.location() << endl
					 << "accepted count: " << (cnt - 1) << endl;

			auto scores = CalculateScoresForWater(water, prot, xmap, distances, bonds, resolution, meanDensity, rmsDensity, *notAHBond, params);
			
			if (dump.is_open())
			{
				dump << cnt << '\t'
					 << id << '\t';
				
				for (double s: scores)
					dump << '\t' << s;

				dump << endl;
			}
			
			if (includeDisputed)
			{
				for (size_t i = 1; i <= proteins.size(); ++i)
				{
					labels.push_back(i <= cnt ? 1 : 0);
					data.push_back(svm::Vector(scores.begin(), scores.end()));
				}
			}
			else if (cnt == 1 or cnt == proteins.size()) 
			{
				labels.push_back(cnt == proteins.size() ? 1 : 0);
				data.push_back(svm::Vector(scores.begin(), scores.end()));
			}
			
			progress.consumed(1);
		}
	}
	
	// print out statistics
	
	cout << "ID" << '\t'
		 << "Water" << '\t'
		 << "Not Water" << '\t'
		 << "Disputed" << '\t'
		 << "Accepted by" << endl;
	for (auto& st: stats)
	{
		cout << st.structure << '\t'
			 << st.accepted << '\t'
			 << st.rejected << '\t'
			 << st.disputed;
		
		for (auto& u: st.acceptedBy)
			cout << '\t' << u.first << ':' << u.second;

		cout << endl;
	}

	auto scales = data.Scale();
	
	if (scales.empty())
		throw runtime_error("Nothing learned");

	if (VERBOSE)
	{
		cout << string(cif::get_terminal_width(), '=') << endl
			 << "Number of waters processed: " << data.size() << endl
			 << endl
			 << "Scales: " << endl;
		
		for (auto& scale: scales)
			cout << '(' << scale.first << ", " << scale.second << ')' << endl;
	}

	typedef svm::SVM_C_SVC_RBF SVM;
	typedef SVM::param_type SVMParams;

	SVMParams svmParams;
	svmParams.probability = true;

	// grid search for C and gamma 
	if (vm.count("grid-search"))
	{
		const int cBegin = -5, cEnd = 11, cStep = 2, gBegin = 3, gEnd = -9, gStep = -2;
		
		vector<pair<SVMParams,size_t>> paramScores;
		for (auto c = cBegin; c <= cEnd; c+= cStep)
			for (auto g = gBegin; g >= gEnd; g += gStep)
			{
				SVMParams svmParams;
				svmParams.gamma = pow(2.0, g);
				svmParams.C = pow(2.0, c);
				
				paramScores.push_back(make_pair(svmParams, 0));
			} 
		
		// use as much threads as possible
		cif::Progress progress(paramScores.size(), "SVM params calculation");
		
		boost::thread_group t;
		size_t N = boost::thread::hardware_concurrency();
		atomic<size_t> next(0);
	
		float bestScoreSoFar = 0;
	
		for (size_t i = 0; i < N; ++i)
			t.create_thread([&]()
			{
				for (;;)
				{
					size_t i = next++;
					
					if (i >= paramScores.size())
						break;

					stringstream s;
					s << "C:" << paramScores[i].first.C << " g:" << paramScores[i].first.gamma;
					progress.message(s.str());
					
					SVM svm(paramScores[i].first);
					auto cv = svm.CrossValidation(labels, data, 5);
					
					size_t tp = 0, tn = 0, fp = 0, fn = 0;

					for (size_t i = 0; i < labels.size(); ++i)
					{
						if (cv[i] == labels[i])
						{
							if (cv[i] == 1)
								++tp;
							else
								++tn;
						}
						else
						{
							if (cv[i] == 1)
								++fp;
							else
								++fn;
						}
					}
		
					paramScores[i].second = tp + tn;

					auto mcc = ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
					
					if (bestScoreSoFar < mcc)
					{
						bestScoreSoFar = mcc;
						cout << endl
							 << setprecision(3) << bestScoreSoFar << " g " << paramScores[i].first.gamma << " C " << paramScores[i].first.C
							 << endl;
					}
	
					progress.consumed(1);
				}
			});
		
		t.join_all();
	
		sort(paramScores.begin(), paramScores.end(),[](const pair<SVMParams,size_t>& a, const pair<SVMParams,size_t>& b) -> bool
			{ return a.second > b.second; });
			
		cout << "Best C: " << paramScores.front().first.C << endl
			 << "Best gamma: " << paramScores.front().first.gamma << endl;
		
		svmParams = paramScores.front().first;
	}
	else
	{
		svmParams.C = 32;
		svmParams.gamma = 8;
	}

	svmParams.probability = true;

	SVM svm(svmParams);

	auto cv = svm.CrossValidation(labels, data, 5);
	
	size_t tp = 0, tn = 0, fp = 0, fn = 0;

	for (size_t i = 0; i < labels.size(); ++i)
	{
		if (cv[i] == labels[i])
		{
			if (cv[i] == 1)
				++tp;
			else
				++tn;
		}
		else
		{
			if (cv[i] == 1)
				++fp;
			else
				++fn;
		}
	}

	auto mcc = ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
	float score = 100.0 * (tp + tn) / labels.size();

	cout << "Cross validation Accuracy: " << score << '%' << endl
		 << "Matthews correlation coefficient: " << mcc << endl
		 << "tp: " << tp << " tn: " << tn << " fp: " << fp << " fn: " << fn << endl;

	auto model = svm.Train(labels, data);
	
	// Store the model and accompanying data
	zeep::xml::document doc("<centrifuge-config/>");
	
	auto config = new zeep::xml::element("params");
	zeep::xml::serializer src(config);
	params.serialize(src, 0);
	doc.child()->append(config);
	
	for (auto& scale: scales)
	{
		auto e = new zeep::xml::element("scale");
		e->set_attribute("min", to_string(scale.first));
		e->set_attribute("max", to_string(scale.second));
		doc.child()->append(e);
	}
	
	doc.child()->append(model->GetConfig());
	
	notAHBond->Save(*doc.child());
	
	fs::ofstream modelFile(vm["output"].as<string>());
	if (not modelFile.is_open())
		throw runtime_error("Could not write to model file");
	
	zeep::xml::writer w(modelFile, true);
	doc.write(w);
	
	return 0;
}

// --------------------------------------------------------------------

int centrifugePredict(int argc, char* argv[])
{
	po::options_description desc("centrifuge-predict " + VERSION + " options");
	desc.add_options()
		("model,m",		po::value<string>(),	"Model file")
		("input-coordinates,i",
						po::value<string>(),	"Input file")
		("input-map,M",	po::value<string>(),	"Input Map or MTZ file")
		("output,o",	po::value<string>(),	"Output file")
		("mtz-columns",
						po::value<string>(),	"Comma separated list of MTZ Column labels to load (default value is taken from model file)")
		("coot-script",	po::value<string>(),	"Write a coot script to this file")
//		("report,r",	po::value<string>(),	"Write the statics report to this file, use 'stdout' to output to screen")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("sampling-rate",
						po::value<float>(),		"Sampling rate")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::positional_options_description p;
	p.add("model", 1);
	p.add("input-coordinates", 1);
	p.add("input-map", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("model") == 0 or vm.count("input-coordinates") == 0)
	{
		cerr << desc << endl;
		exit(1);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();
	
	// Read the model
	fs::ifstream modelFile(vm["model"].as<string>());
	if (not modelFile.is_open())
		throw runtime_error("no such model file");
	
	zeep::xml::document doc(modelFile);
	auto e = doc.child()->find_first("params");
	if (e == nullptr)
		throw runtime_error("invalid model file");
	
	CentrifugeParameters params;
	zeep::xml::deserializer dsp(e);
	params.serialize(dsp, 0);
	
	vector<pair<float,float>> scales;
	for (auto e: doc.child()->find("scale"))
	{
		scales.push_back(make_pair(
			stof(e->get_attribute("min")),
			stof(e->get_attribute("max"))));
	}
	
	if (scales.size() != kNrOfFeatures)
		throw runtime_error("Invalid number of scales in model file");

	e = doc.child()->find_first("svm-config");
	if (e == nullptr)
		throw runtime_error("invalid model file");
	auto model = svm::ModelBase::Create(*e);
	
	fs::path file = vm["input-coordinates"].as<string>();

	if (VERBOSE)
		cout << string(cif::get_terminal_width(), '=') << endl
			 << "Processing file: " << file << endl;
	
	fs::ofstream cootScript;
	if (vm.count("coot-script"))
	{
		cootScript.open(vm["coot-script"].as<string>());
		if (not cootScript.is_open())
			throw runtime_error("Could not create coot script");
		
		cootScript << "(interesting-things-gui \"Waters\"" << endl
					<< "(list" << endl;
	}
	
	c::File mmcif(file);
	c::Structure prot(mmcif, 1);

//	clipper::Xmap<float> xmap;
//	float resolution = db["refine"].front()["ls_d_res_high"].as<float>(), meanDensity, rmsDensity;
//	fs::path mapFile;
//
//	auto mtzColumns = params.mtzColumns;
//	if (vm.count("mtz-columns"))
//	{
//		mtzColumns.clear();
//		ba::split(mtzColumns, vm["mtz-columns"].as<string>(), ba::is_any_of(","));
//	}
//
//	if (vm.count("input-map"))
//		mapFile = vm["input-map"].as<string>();
//	else
//	{
//		if (file.extension() == ".gz" or file.extension() == ".bz2")
//			file = file.parent_path() / file.stem();
//		
//		mapFile = file.parent_path() / (file.stem().string() + ".ccp4");
//		if (not fs::exists(mapFile))
//			mapFile = file.parent_path() / (file.stem().string() + ".ccp4.gz");
//		if (not fs::exists(mapFile))
//			mapFile = file.parent_path() / (file.stem().string() + ".ccp4.bz2");
//		if (not fs::exists(mapFile))
//			mapFile = file.parent_path() / (file.stem().string() + ".mtz");
//		if (not fs::exists(mapFile))
//			mapFile = file.parent_path() / (file.stem().string() + ".mtz.gz");
//		if (not fs::exists(mapFile))
//			mapFile = file.parent_path() / (file.stem().string() + ".mtz.bz2");
//	}
//
//	ReadXmap(mapFile, mtzColumns, xmap, resolution, meanDensity, rmsDensity);

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

	c::MapMaker<float> mm;
	mm.loadMTZ(vm["input-map"].as<string>(), samplingRate);

	clipper::Xmap<float>& xmap = mm.fb();
	
	float resolution = mm.resHigh(), meanDensity = mm.fb().meanDensity(), rmsDensity = mm.fb().rmsDensity();

	// Calculate the distance map
	auto& spacegroup = xmap.spacegroup();
	auto& cell = xmap.cell();
	
	c::DistanceMap distances(prot, spacegroup, cell, 3.5f);
	c::BondMap bonds(prot);
	
	// Load the NotAHBondSet
	auto notAHBond = NotAHBondSet::Load(*doc.child());
	vector<c::Atom> notWater;
	
	for (auto& water : prot.waters())
	{
		auto scores = CalculateScoresForWater(water, prot, xmap, distances, bonds, resolution, meanDensity, rmsDensity, *notAHBond, params);
		
		svm::Vector v;
		
		for (size_t i = 0; i < kNrOfFeatures; ++i)
			v[i] = (scores[i] - scales[i].first) / (scales[i].second - scales[i].first);

		// score with probability			
		auto p = model->PredictWithProbability(v);
		
		// score without probability
		bool isWater = model->Predict(v);
		
		auto voted = distance(p.begin(), max_element(p.begin(), p.end()));
		
		cout << water.id() << '\t'
			 << water.authAsymId() << ':' << water.authSeqId() << "\t"
			 << (isWater ? "water" : "not water") << '\t'
			 << '(' << (p[0] * 100) << "%)"
			 << (voted == isWater ? "!" : "")
			 << endl;
		
		if (cootScript.is_open())
		{
			float x, y, z;
			tie(x, y, z) = water.location();
			cootScript << "(list \"" << (isWater ? "Keep" : "Remove")
						<< " HOH " << water.id()
						<< " (" << setprecision(3) <<  (p[0] * 100) << "%)\" "
						<< x << ' ' << y << ' ' << z << " )" << endl;
		}
		
		if (not isWater)
			 notWater.push_back(water);
	}
	
	if (vm.count("output"))
	{
		for (auto& a: notWater)
			prot.removeAtom(a);
		
		prot.getFile().save(vm["output"].as<string>());
	}
	
	if (cootScript.is_open())
		cootScript << "))" << endl;

	return 0;
}

// --------------------------------------------------------------------

int edia_test(int argc, char* argv[])
{
	po::options_description desc("centrifuge-predict " + VERSION + " options");
	desc.add_options()
		("pdb",			po::value<string>(),	"Input PDB file")
		("mtz",			po::value<string>(),	"Input MTZ or CCP4MAP file")
		("output,o",	po::value<string>(),	"Output file")
		("mtz-columns",
						po::value<string>(),	"Comma separated list of MTZ Column labels to load (default value is \"FWT,PHWT\"")
		("resolution",	po::value<float>(),		"Resolution to use, default is taken from mtz file or pdb file")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		("sampling-rate",
						po::value<float>(),		"Sampling rate")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::positional_options_description p;
	p.add("pdb", 1);
	p.add("mtz", 2);
	p.add("output", 3);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("pdb") == 0 or vm.count("mtz") == 0)
	{
		cerr << desc << endl;
		exit(1);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();

	vector<string> mtzColumns = { "FWT", "PHWT" };

	if (vm.count("mtz-columns"))
		ba::split(mtzColumns, vm["mtz-columns"].as<string>(), ba::is_any_of(","));
	
	fs::path pdbFile = vm["pdb"].as<string>();
	fs::path mtzFile = vm["mtz"].as<string>();

	if (VERBOSE)
		cout << string(cif::get_terminal_width(), '=') << endl
			 << "Processing file: " << pdbFile << endl;
	
	c::File mmcif(pdbFile);
	c::Structure prot(mmcif, 1);
	cif::Datablock& db = mmcif.data();

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

	c::MapMaker<float> mm;
	mm.loadMTZ(mtzFile, samplingRate);

	clipper::Xmap<float>& xmap = mm.fb();
	
	float resolution = mm.resHigh(), meanDensity = mm.fb().meanDensity(), rmsDensity = mm.fb().rmsDensity();

	if (vm.count("resolution"))
		resolution = vm["resolution"].as<float>();

	if (VERBOSE)
		cout << "Resolution: " << resolution << endl
			 << "Mean Density: " << meanDensity << endl
			 << "RMS Density: " << rmsDensity << endl
			 << endl;

	// Calculate the distance map
	auto& spacegroup = xmap.spacegroup();
	auto& cell = xmap.cell();
	
	c::DistanceMap distances(prot, spacegroup, cell, 3.5f);
	c::BondMap bonds(prot);

	map<string,float> ediaScores;

//	for (auto& water : prot.waters())
	for (auto& a : prot.atoms())
	{
//		float edia = CalculateEDIA(a, xmap, resolution, meanDensity, rmsDensity, distances, bonds);
//		
//		ediaScores[a.id()] = edia;
//		
//		cout << a.authAsymId() << ':' << a.authSeqId() << '\t'
//			 << a.labelAtomId() << '\t'
//			 << c::AtomTypeTraits(a.type()).symbol() << '\t'
//			 << fixed << setprecision(2) << edia << endl;
	}
	
	// And EDIAm scores

	cout << endl
		 << string(cif::get_terminal_width(), '=') << endl
		 << "EDIAm scores for resdidues: " << endl
		 << "Asym/seqID" << '\t'
		 << "comp" << '\t'
		 << "EDIAm" << '\t'
		 << "min EDIA" << '\t'
		 << "median EDIA" << '\t'
		 << "max EDIA" << '\t'
		 << "OPIA" << endl;

	vector<tuple<string,string,int>> residues;

	const int kNonPolySeqID = numeric_limits<int>::max();

	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		string asymId, compId;
		int seqId;
		
		cif::tie(asymId, compId, seqId) = r.get("asym_id", "mon_id", "seq_id");
		residues.push_back(make_tuple(asymId, compId, seqId));
	}

	for (auto r: db["pdbx_nonpoly_scheme"])
	{
		string asymId, compId;
		
		cif::tie(asymId, compId) = r.get("asym_id", "mon_id");
		residues.push_back(make_tuple(asymId, compId, kNonPolySeqID));
	}
	
	for (size_t i = 0; i < residues.size(); ++i)
	{
		string asymId, compId;
		int seqId;
		
		tie(asymId, compId, seqId) = residues[i];
		
		auto compound = mmcif::Compound::create(compId);
		if (not compound or compound->isWater())
			continue;
		
		vector<double> edia;
		
		double sum = 0;
		size_t N = 0, M = 0;
		
		vector<string> compoundAtoms;
		if (seqId != kNonPolySeqID and (i + 1 == residues.size() or get<0>(residues[i + 1]) != asymId))
			compoundAtoms = { "OXT" };

		for (auto cAtom: compound->atoms())
		{
			if (cAtom.typeSymbol == mmcif::H)
				continue;
			compoundAtoms.push_back(cAtom.id);
		}
		
		for (auto cAtomID: compoundAtoms)
		{
			++N;
			
			auto a = 
				seqId != kNonPolySeqID ?
					db["atom_site"].find(cif::Key("label_asym_id") == asymId and cif::Key("label_seq_id") == seqId and cif::Key("label_atom_id") == cAtomID) :
					db["atom_site"].find(cif::Key("label_asym_id") == asymId and cif::Key("label_atom_id") == cAtomID);

			if (a.size() > 1)
				throw runtime_error("Invalid number of atoms found for asym_id " + asymId + " comp_id " + compId + (seqId == kNonPolySeqID ? "" : " seq_id " + to_string(seqId)));

			if (a.empty())
			{
				if (cAtomID == "OXT")
					--N;
				else if (VERBOSE)
					cerr << "Missing atom for compound " << compId << " in asym " << asymId << " seq_id " << (seqId == kNonPolySeqID ? 0 : seqId) << endl;
			}
			else
			{
				string atomId = a.front()["id"].as<string>();
				
				float EDIAa = ediaScores[atomId];
				
				if (EDIAa >= 0.8)
					++M;

				edia.push_back(EDIAa);
				sum += pow(EDIAa + 0.1, -2);
			}
		}
		
		float EDIAm = 0, ediaMin = 0, ediaMax = 0, ediaMedian = 0, OPIA = 0;
		
		if (not edia.empty())
		{
			EDIAm = static_cast<float>(pow(sum / N, -0.5) - 0.1);
			
			sort(edia.begin(), edia.end());
	
			ediaMin = edia.front();
			ediaMax = edia.back();
			
			if (edia.size() % 1)
				ediaMedian = edia[edia.size() / 2];
			else
				ediaMedian = (edia[edia.size() / 2] + edia[edia.size() / 2 + 1]) / 2;
		
			OPIA = 100.0f * M / N;
		}
		
		cout << asymId << '\t'
			 << compId << '\t'
			 << setprecision(2) << EDIAm << '\t'
			 << setprecision(2) << ediaMin << '\t'
			 << setprecision(2) << ediaMedian << '\t'
			 << setprecision(2) << ediaMax << '\t'
			 << setprecision(1) << OPIA << endl;
	}
	
	return 0;
}

// --------------------------------------------------------------------
//    main

int pr_main(int argc, char* argv[])
{
	int result = 0;
	
	string command = fs::path(argv[0]).leaf().string();

	if (command == "centrifuge" and argc > 1)
	{
		command = argv[1];
		++argv;
		--argc;
	}
	
	if (command == "learn")
		result = centrifugeLearn(argc, argv);
	else if (command == "predict")
		result = centrifugePredict(argc, argv);

// other commands here
	else
	{
		cout << "Usage: centrifuge command [options]" << endl
			 << endl
			 << "Where command is either learn or predict" << endl
			 << endl;

		exit(1);
	}
	
	return result;
}
