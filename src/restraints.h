/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 22 mei, 2018
*/

#pragma once

#include <set>

#include "cif++/Point.hpp"
#include "cif++/MapMaker.hpp"

// --------------------------------------------------------------------

class AtomLocationProvider;
class DFCollector;

// --------------------------------------------------------------------

typedef uint32_t AtomRef;
typedef typename mmcif::Map<float>::Xmap Xmap;

// --------------------------------------------------------------------

struct Restraint
{
	virtual ~Restraint() {}
	
	virtual double f(const AtomLocationProvider& atoms) const = 0;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const = 0;
	
	virtual void print(const AtomLocationProvider& atoms) const = 0;
};

struct BondRestraint : public Restraint
{
	BondRestraint(AtomRef a, AtomRef b, double distance, double esd)
		: mA(a), mB(b), mDist(distance), mDistESD(esd) {}

	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;
	
	AtomRef	mA, mB;
	double	mDist, mDistESD;
};

struct AngleRestraint : public Restraint
{
	AngleRestraint(AtomRef a, AtomRef b, AtomRef c, double angle, double esd)
		: mA(a), mB(b), mC(c), mAngle(angle), mESD(esd) {}
	
	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;

	AtomRef	mA, mB, mC;
	double	mAngle, mESD;
};

struct TorsionRestraint : public Restraint
{
	TorsionRestraint(AtomRef a, AtomRef b, AtomRef c, AtomRef d, double target, double esd, int periodicity)
		: mA(a), mB(b), mC(c), mD(d), mPeriodicity(periodicity), mTarget(target), mESD(esd) {}
	
	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;

	AtomRef	mA, mB, mC, mD;
	int mPeriodicity;
	double	mTarget, mESD;

  private:
	std::tuple<mmcif::DPoint,mmcif::DPoint,mmcif::DPoint,mmcif::DPoint> CalculateTorsionGradients(float theta, mmcif::DPoint p[4]) const;
};

struct TransPeptideRestraint : public TorsionRestraint
{
	TransPeptideRestraint(AtomRef a, AtomRef b, AtomRef c, AtomRef d, double esd = 2.0)
		: TorsionRestraint(a, b, c, d, 180.0, esd, 2) {}
};

struct ChiralVolumeRestraint : public Restraint
{
	ChiralVolumeRestraint(AtomRef c, AtomRef a1, AtomRef a2, AtomRef a3, double volume)
		: mCentre(c), mA1(a1), mA2(a2), mA3(a3), mVolume(volume) {}
	
	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;

	AtomRef	mCentre, mA1, mA2, mA3;
	double	mVolume;
};

struct PlanarityRestraint : public Restraint
{
	PlanarityRestraint(std::vector<AtomRef>&& atoms, double esd)
		: mAtoms(std::move(atoms)), mESD(esd)
	{
		if (mAtoms.size() < 3)
			throw std::runtime_error("Insufficient number of atoms in planar restraint");
	}
	
	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;

	void calculatePlaneFunction(const AtomLocationProvider& atoms, double abcd[4]) const;
	
	std::vector<AtomRef> mAtoms;
	double mESD;
};

struct NonBondedContactRestraint : public Restraint
{
	NonBondedContactRestraint(AtomRef a, AtomRef b, double minDist, double esd)	
		: mA(a), mB(b), mMinDist(minDist)
		, mMinDistSq(minDist * minDist), mDistESD(esd) {}
	
	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;

	AtomRef	mA, mB;
	double	mMinDist, mMinDistSq, mDistESD;
};

struct DensityRestraint : public Restraint
{
	DensityRestraint(std::vector<std::pair<AtomRef,double>>&& atoms,
		const Xmap& xMap, double mapWeight = 60);

	virtual double f(const AtomLocationProvider& atoms) const;
	virtual void df(const AtomLocationProvider& atoms, DFCollector& d) const;
	virtual void print(const AtomLocationProvider& atoms) const;
	
	std::vector<std::pair<AtomRef,double>> mAtoms;
	const Xmap& mXMap;
	double mMapWeight;
	bool mElectronScattering = false;
};
