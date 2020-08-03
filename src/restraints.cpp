/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 22 mei, 2018
*/

#include "pdb-redo.h"

#include <numeric>

#include <boost/algorithm/string.hpp>

#include "cif++/Structure.hpp"
#include "cif++/CifUtils.hpp"

#include "restraints.h"
#include "minimizer.h"

using namespace std;

namespace ba = boost::algorithm;

using mmcif::Atom;
using mmcif::Point;
using mmcif::DPoint;
using mmcif::Structure;
using mmcif::Monomer;
using mmcif::Distance;

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Coord_map;
using clipper::Coord_frac;

// --------------------------------------------------------------------

double BondRestraint::f(const AtomLocationProvider& atoms) const
{
	double d = mDist - Distance(atoms[mA], atoms[mB]);
	double result = (d * d) / (mDistESD * mDistESD);
	
	if (cif::VERBOSE > 2)
		cerr << "bond::f() = " << atoms.atom(mA) << " <> " << atoms.atom(mB)
			 << " " << atoms[mA] << " <> " << atoms[mB]
			 << " => " << result << endl;
	
	return result;
}

void BondRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	auto a1 = atoms[mA], a2 = atoms[mB];
	
	auto bi = Distance(a1, a2);
	if (bi < 0.1)
		bi = 0.1;
	
	auto c = 2 * (1 - mDist / bi) / (mDistESD * mDistESD);
	
	if (cif::VERBOSE > 2)
		cerr << "bond::df(): " << atoms.atom(mA) << " <> " << atoms.atom(mB) << ' '
			 << bi << ' ' << mDist << ' ' << mDistESD << endl;

	df.add(mA, (a1 - a2) * c);
	df.add(mB, (a2 - a1) * c);
}

void BondRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "bond " << atoms.atom(mA) << " to " << atoms.atom(mB) << " => " << mDist << " / " << mDistESD << endl;
}

// --------------------------------------------------------------------

double AngleRestraint::f(const AtomLocationProvider& atoms) const
{
	DPoint p[3] = { atoms[mA], atoms[mB], atoms[mC] };
	
	double c = CosinusAngle(p[1], p[0], p[1], p[2]);
	double angle = atan2(sqrt(1 - c * c), c) * 180 / mmcif::kPI;

	double d = mAngle - angle;
	double result = (d * d) / (mESD * mESD);
	
	if (cif::VERBOSE > 2)
		cerr << "angle::f() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << ' ' << " = " << result << endl;
	
	return result;
}

void AngleRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	const double kRadToDegree = 180.0 / mmcif::kPI, kDegreeToRad = 1 / kRadToDegree;
	
	if (cif::VERBOSE > 2)
		cerr << "angle::df() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << ' ' << ": " << endl;

	DPoint k = atoms[mA], l = atoms[mB], m = atoms[mC];
	
	auto aVec = (k - l);
	auto bVec = (m - l);
	
	auto a = Distance(k, l);
	if (a < 0.01)
	{
		a = 0.01;
		aVec = DPoint{ 0.01, 0.01, 0.01 };
	};
		
	auto b = Distance(m, l);
	if (b < 0.01)
	{
		b = 0.01;
		bVec = DPoint{ 0.01, 0.01, 0.01 };
	};
	
	auto cosTheta = DotProduct(aVec, bVec) / (a * b);
	if (cosTheta > 1.0)	cosTheta = 1.0;
	if (cosTheta < -1.0) cosTheta = -1.0;
	auto theta = acos(cosTheta);
	if (theta < 0.001) theta = 0.001;
	
	auto target = mAngle * kDegreeToRad;
	
	auto wf = 2 * (theta - target) * kRadToDegree * kRadToDegree / (mESD * mESD);
	auto prem = -wf / sin(theta);
	
	df.add(mA, prem * (cosTheta * (l - k) / (a * a) + (m - l) / (a * b)));
	df.add(mC, prem * (cosTheta * (l - m) / (b * b) + (k - l) / (a * b)));
	
	auto term1 = (l - k) * -cosTheta / (a * a) + (l - m) * -cosTheta / (b * b);
	auto term2 = ((l - k) + (l - m)) / (a * b);
	
	df.add(mB, prem * (term1 + term2));
}

void AngleRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "angle " << atoms.atom(mA) << " ; " << atoms.atom(mB) << " ; " << atoms.atom(mC) << " => " << mAngle << " / " << mESD << endl;
}

// --------------------------------------------------------------------

tuple<DPoint,DPoint,DPoint,DPoint>
TorsionRestraint::CalculateTorsionGradients(float theta, DPoint p[4]) const
{
	auto a = p[1] - p[0], b = p[2] - p[1], c = p[3] - p[2];
	
	auto blensq = b.lengthsq();
	auto blen = sqrt(blensq);
	
	if (blen < 0.01)
	{
		blen = 0.01;
		blensq = 0.0001;
	}
	
	auto H = -DotProduct(a, c),
		 J =  DotProduct(a, b),
		 K =  DotProduct(b, c),
		 L = 1 / blensq;
	
	auto E = DotProduct(a, CrossProduct(b, c)) / blen;
	auto G = H + J * K * L;
	auto F = 1 / G;
	
	if (G == 0)
		F = 999999999.9;
	
	DPoint dH[4] = { c, -c, a, -a };
	DPoint dK[4] = { {}, -c, c - b, b };
	DPoint dJ[4] = { -b, b - a, a, {} };
	DPoint dL[4] = { {}, 2.0 * (p[2] - p[1]) * L * L, -2.0 * (p[2] - p[1]) * L * L, {} };
	DPoint dM[4] = {
		{
		   -(b.mY * c.mZ - b.mZ * c.mY),
		   -(b.mZ * c.mX - b.mX * c.mZ),
		   -(b.mX * c.mY - b.mY * c.mX)
		},
		{
			(b.mY * c.mZ - b.mZ * c.mY) + (a.mY * c.mZ - a.mZ * c.mY),
			(b.mZ * c.mX - b.mX * c.mZ) + (a.mZ * c.mX - a.mX * c.mZ),
			(b.mX * c.mY - b.mY * c.mX) + (a.mX * c.mY - a.mY * c.mX)
		},
		{
			(b.mY * a.mZ - b.mZ * a.mY) - (a.mY * c.mZ - a.mZ * c.mY),
		   -(a.mZ * c.mX - a.mX * c.mZ) + (b.mZ * a.mX - b.mX * a.mZ),
		   -(a.mX * c.mY - a.mY * c.mX) + (a.mY * b.mX - a.mX * b.mY)
		},
		{
		   -(b.mY * a.mZ - b.mZ * a.mY),
		   -(b.mZ * a.mX - b.mX * a.mZ),
		   -(a.mY * b.mX - a.mX * b.mY)
		}
	};
	
	DPoint dE[4] {
		dM[0] / blen,
		dM[1] / blen + E * (p[2] - p[1]) * L,
		dM[2] / blen - E * (p[2] - p[1]) * L,
		dM[3] / blen
	};
	
	auto eff = E * F * F;
	auto jl = J * L;
	auto kl = K * L;
	auto jk = J * K;
	
	return make_tuple(
		F * dE[0] - eff * (dH[0] + jl * dK[0] + kl * dJ[0] + jk * dL[0]),
		F * dE[1] - eff * (dH[1] + jl * dK[1] + kl * dJ[1] + jk * dL[1]),
		F * dE[2] - eff * (dH[2] + jl * dK[2] + kl * dJ[2] + jk * dL[2]),
		F * dE[3] - eff * (dH[3] + jl * dK[3] + kl * dJ[3] + jk * dL[3])
	);
}

double TorsionRestraint::f(const AtomLocationProvider& atoms) const
{
	double result = 0;
	
	double cos_a1 = CosinusAngle(atoms[mB], atoms[mA], atoms[mC], atoms[mB]);
	double cos_a2 = CosinusAngle(atoms[mC], atoms[mB], atoms[mD], atoms[mC]);
	
	if (cos_a1 <= 0.9 and cos_a2 <= 0.9)
	{
		double period = 360.0;
		if (mPeriodicity > 0)
			period /= mPeriodicity;
		
		double theta = DihedralAngle(atoms[mA], atoms[mB], atoms[mC], atoms[mD]);
		double diff = fmod(abs(theta - mTarget) + period / 2, period) - period / 2;
	
		if (not isnan(diff))
			result = (diff * diff) / (mESD * mESD);

		if (cif::VERBOSE > 2)
			cerr << "torsion::f() = " << result << " for theta " << theta
				 << " diff: " << diff
				 << " target: " << mTarget << " sigma: " << mESD
				 << " atoms: " << atoms.atom(mA) << ", " << atoms.atom(mB) << ", " << atoms.atom(mC) << ", " << atoms.atom(mD)
				 << endl;
	}
	
	return result;
}

void TorsionRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	if (cif::VERBOSE > 2)
		cerr << "torsion::df() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << "/" << atoms.atom(mD) << ' ' << ": " << endl;

	double cos_a1 = CosinusAngle(atoms[mB], atoms[mA], atoms[mC], atoms[mB]);
	double cos_a2 = CosinusAngle(atoms[mC], atoms[mB], atoms[mD], atoms[mC]);
	
	if (cos_a1 <= 0.9 and cos_a2 <= 0.9)
	{
		double period = 360.0;
		if (mPeriodicity > 0)
			period /= mPeriodicity;
		
		double theta = DihedralAngle(atoms[mA], atoms[mB], atoms[mC], atoms[mD]);
		double diff = theta - mTarget;
		if (diff > 180)
			diff -= 360;
		if (diff < -180)
			diff += 360;
	
		if (not isnan(diff))
		{
			auto tt = tan(mmcif::kPI * theta / 180);
			double scale = 180.0 / ((1 + tt * tt) * mmcif::kPI);
			auto w = 1 / (mESD * mESD);
			
			DPoint p[4] = { atoms[mA], atoms[mB], atoms[mC], atoms[mD] };
			DPoint d[4];
			
			tie(d[0], d[1], d[2], d[3]) = CalculateTorsionGradients(theta, p);

			df.add(mA, 2.0 * diff * d[0] * scale * w);
			df.add(mB, 2.0 * diff * d[1] * scale * w);
			df.add(mC, 2.0 * diff * d[2] * scale * w);
			df.add(mD, 2.0 * diff * d[3] * scale * w);
		}
	}	
}

void TorsionRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "torsion " << atoms.atom(mA) << " ; " << atoms.atom(mB) << " ; " << atoms.atom(mC) << " ; " << atoms.atom(mD) << " => " << mPeriodicity << " / " << mESD << endl;
}

// --------------------------------------------------------------------

const double kChiralVolumeESD = 0.2;	// according to coot that's a reasonable value...

double ChiralVolumeRestraint::f(const AtomLocationProvider& atoms) const
{
	auto chiralVolume = DotProduct(atoms[mA1] - atoms[mCentre],
		CrossProduct(atoms[mA2] - atoms[mCentre], atoms[mA3] - atoms[mCentre]));

	double d = mVolume - chiralVolume;
	double result = (d * d) / (kChiralVolumeESD * kChiralVolumeESD);
	
	if (cif::VERBOSE > 2)
		cerr << "chiral::f() = " << result << endl;
	
	return result;
}

void ChiralVolumeRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	if (cif::VERBOSE > 2)
		cerr << "chiral::df(): " << endl;

	DPoint centre = atoms[mCentre];
	DPoint a = atoms[mA1] - centre;
	DPoint b = atoms[mA2] - centre;
	DPoint c = atoms[mA3] - centre;
	
	auto chiralVolume = DotProduct(a, CrossProduct(b, c));
	
	auto d = chiralVolume - mVolume;
	auto s = 2 * d / (kChiralVolumeESD * kChiralVolumeESD);
	
	df.add(mCentre, s * DPoint{
		- (b.mY * c.mZ - b.mZ * c.mY) - (a.mZ * c.mY - a.mY * c.mZ) - (a.mY * b.mZ - a.mZ * b.mY),
		- (b.mZ * c.mX - b.mX * c.mZ) - (a.mX * c.mZ - a.mZ * c.mX) - (a.mZ * b.mX - a.mX * b.mZ),
		- (b.mX * c.mY - b.mY * c.mX) - (a.mY * c.mX - a.mX * c.mY) - (a.mX * b.mY - a.mY * b.mX)
	});

	df.add(mA1, s * DPoint{ b.mY * c.mZ - b.mZ * c.mY, b.mZ * c.mX - b.mX * c.mZ, b.mX * c.mY - b.mY * c.mX});
	df.add(mA2, s * DPoint{ a.mZ * c.mY - a.mY * c.mZ, a.mX * c.mZ - a.mZ * c.mX, a.mY * c.mX - a.mX * c.mY});
	df.add(mA3, s * DPoint{ a.mY * b.mZ - a.mZ * b.mY, a.mZ * b.mX - a.mX * b.mZ, a.mX * b.mY - a.mY * b.mX});
}

void ChiralVolumeRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "chiral volume " << atoms.atom(mA1) << " ; " << atoms.atom(mA2) << " ; " << atoms.atom(mA3) << " => " << mVolume << " / " << kChiralVolumeESD << endl;
}

// --------------------------------------------------------------------

void PlanarityRestraint::calculatePlaneFunction(const AtomLocationProvider& atoms, double abcd[4]) const
{
	DPoint center;
	
	for (auto& a: mAtoms)
		center += atoms[a];
	center /= mAtoms.size();
	
	clipper::Matrix<double> mat(3, 3);
	for (auto& a: mAtoms)
	{
		mat(0, 0) += (atoms[a].mX - center.mX) * (atoms[a].mX - center.mX);
		mat(1, 1) += (atoms[a].mY - center.mY) * (atoms[a].mY - center.mY);
		mat(2, 2) += (atoms[a].mZ - center.mZ) * (atoms[a].mZ - center.mZ);
		mat(0, 1) += (atoms[a].mX - center.mX) * (atoms[a].mY - center.mY);
		mat(0, 2) += (atoms[a].mX - center.mX) * (atoms[a].mZ - center.mZ);
		mat(1, 2) += (atoms[a].mY - center.mY) * (atoms[a].mZ - center.mZ);
	}
	
	mat(1, 0) = mat(0, 1);
	mat(2, 0) = mat(0, 2);
	mat(2, 1) = mat(1, 2);
	
	auto eigen = mat.eigen(true);
	
	abcd[0] = mat(0, 0);
	abcd[1] = mat(1, 0);
	abcd[2] = mat(2, 0);

	double sumSq = 1e-20 + abcd[0] * abcd[0] + abcd[1] * abcd[1] + abcd[2] * abcd[2];

	abcd[0] /= sumSq;
	abcd[1] /= sumSq;
	abcd[2] /= sumSq;
	
	abcd[3] = abcd[0] * center.mX + abcd[1] * center.mY + abcd[2] * center.mZ;
}

double PlanarityRestraint::f(const AtomLocationProvider& atoms) const
{
	double abcd[4];

	calculatePlaneFunction(atoms, abcd);
	
	double result = accumulate(mAtoms.begin(), mAtoms.end(), 0.,
		[&atoms,&abcd,this](double sum, AtomRef a)
		{
			double v =	abcd[0] * atoms[a].mX +
						abcd[1] * atoms[a].mY +
						abcd[2] * atoms[a].mZ -
						abcd[3];
			
			double r = v / mESD;
			
			return sum + r * r;
		});

	if (cif::VERBOSE > 2)
	{
		vector<string> as;
		transform(mAtoms.begin(), mAtoms.end(), back_inserter(as),
			[](auto& a)
			{
				stringstream s;
				s << a;
				return s.str();
			});
		
		cerr << "plane::f() = " << result << " for " << mAtoms.size() << " atoms " << ba::join(as, ", ") << endl;
	}
	
	return result;
}

void PlanarityRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	if (cif::VERBOSE > 2)
	{
		vector<string> as;
		transform(mAtoms.begin(), mAtoms.end(), back_inserter(as),
			[](auto& a)
			{
				stringstream s;
				s << a;
				return s.str();
			});
		
		cerr << "plane::df() for " << mAtoms.size() << " atoms " << ba::join(as, ", ") << endl;
	}

	double abcd[4];

	calculatePlaneFunction(atoms, abcd);
	
	for (auto& a: mAtoms)
	{
		auto l = atoms[a];
		auto deviLen = l.mX * abcd[0] + l.mY * abcd[1] + l.mZ * abcd[2] - abcd[3];
		
		df.add(a, 2 * deviLen * DPoint{ abcd[0], abcd[1], abcd[2] } / (mESD * mESD));
	}
}

void PlanarityRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "plane ";
	
	for (auto& a: mAtoms)
		cout << atoms.atom(a) << ' ';

	cout << "=> " << 0 << " / " << mESD << endl;
}

// --------------------------------------------------------------------

double NonBondedContactRestraint::f(const AtomLocationProvider& atoms) const
{
	double result = 0;
	
	double distance = DistanceSquared(atoms[mA], atoms[mB]);
	if (distance < mMinDistSq)
	{
		double d = mMinDist - sqrt(distance);
		result = (d * d) / (mDistESD * mDistESD);
	}
	
	if (cif::VERBOSE > 2)
			cerr << "non-bonded-contact::f() = " << result << " min-dist is " << mMinDist << " and dist is " << sqrt(distance)
			 << " a1: " << atoms.atom(mA) << " a2: " << atoms.atom(mB) << endl
			 << " a1: " << atoms[mA] << " a2: " << atoms[mB] << endl;
	
	return result;
}

void NonBondedContactRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	auto a1 = atoms[mA], a2 = atoms[mB];
	
	auto bi = DistanceSquared(a1, a2);
	if (bi < mMinDistSq)
	{
		bi = sqrt(bi);
		if (bi < 0.1)
			bi = 0.1;
		
		if (cif::VERBOSE > 2)
			cerr << "non-bonded::df(): " << atoms.atom(mA) << " and " << atoms.atom(mB) << " "
				 << "distance: " << bi << " "
				 << "target: " << mMinDist << endl;

		double c = 2 * (1 - mMinDist / bi) / (mDistESD * mDistESD);
		
		df.add(mA, (a1 - a2) * c);
		df.add(mB, (a2 - a1) * c);
	}
}

void NonBondedContactRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "nbc " << atoms.atom(mA) << " " << atoms.atom(mB)
		 << " => " << Distance(atoms[mA], atoms[mB]) << ' ' << mMinDist << " / " << mDistESD << endl;
}

// --------------------------------------------------------------------

DensityRestraint::DensityRestraint(vector<pair<AtomRef,double>>&& atoms,
	const Xmap& xMap, double mapWeight)
	: mAtoms(move(atoms)), mXMap(xMap), mMapWeight(mapWeight)
{
}

double DensityRestraint::f(const AtomLocationProvider& atoms) const
{
	double result = 0;
	
	for (auto& a: mAtoms)
	{
		Coord_orth p = atoms[a.first];
		Coord_frac pf = p.coord_frac(mXMap.cell());
		
		result += a.second * mXMap.interp<clipper::Interp_cubic>(pf);
	}

	if (cif::VERBOSE > 2)
		cerr << "density::f() = " << -result << endl;
	
	return mMapWeight * -result;
}

void DensityRestraint::df(const AtomLocationProvider& atoms, DFCollector& df) const
{
	if (cif::VERBOSE > 2)
		cerr << "density::df(): " << endl;

	for (auto& a: mAtoms)
	{
		Coord_orth p = atoms[a.first];
		Coord_frac pf = p.coord_frac(mXMap.cell());
		auto pm = pf.coord_map(mXMap.grid_sampling());
		
		clipper::Grad_map<double> grad;
		double dv;
		
		clipper::Interp_cubic::interp_grad(mXMap, pm, dv, grad);
		auto gradFrac = grad.grad_frac(mXMap.grid_sampling());
		
		auto gradOrth = gradFrac.grad_orth(mXMap.cell());
		
		df.add(a.first, DPoint{ gradOrth.dx(), gradOrth.dy(), gradOrth.dz() } * mMapWeight * -a.second);
	}
}

void DensityRestraint::print(const AtomLocationProvider& atoms) const
{
	cout << "density " << endl;
}

