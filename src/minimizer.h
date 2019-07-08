/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 22 mei, 2018
*/

#pragma once

#include <boost/optional.hpp>

#include "cif++/Structure.h"
#include "cif++/MapMaker.h"
#include "cif++/BondMap.h"

#include "restraints.h"

// --------------------------------------------------------------------

class AtomLocationProvider
{
  public:
	AtomLocationProvider(const AtomLocationProvider&) = delete;
	AtomLocationProvider& operator=(const AtomLocationProvider&) = delete;

	AtomLocationProvider(std::vector<mmcif::Atom>& atoms)
		: mAtoms(atoms) {}
	virtual ~AtomLocationProvider() {}
	
	virtual mmcif::DPoint operator[](AtomRef atomID) const;
	virtual std::string atom(AtomRef atomID) const;

  protected:
	std::vector<mmcif::Atom>& mAtoms;
};

// --------------------------------------------------------------------

class DFCollector
{
  public:
	DFCollector(const DFCollector&) = delete;
	DFCollector& operator=(const DFCollector&) = delete;

	DFCollector() {}
	virtual ~DFCollector() {}
	
	virtual void add(AtomRef atom, double dx, double dy, double dz) = 0;
	void add(AtomRef atom, mmcif::DPoint&& d)
	{
		add(atom, d.mX, d.mY, d.mZ);
	}
};

// --------------------------------------------------------------------

class Minimizer
{
  public:
	typedef clipper::Xmap<float> XMap;

	Minimizer(const Minimizer&) = delete;
	Minimizer& operator=(const Minimizer&) = delete;

	virtual ~Minimizer() {}

	// factory method:
	static Minimizer* create(const std::string& algorithm,
		const mmcif::Polymer& poly, int first, int last, const mmcif::BondMap& bm,
		const XMap& xMap, float mapWeight = 60, float plane5AtomsESD = 0.11);

	void printStats();

	virtual double refine(bool storeAtoms) = 0;
	double score();
	virtual std::vector<std::pair<std::string,mmcif::Point>> getAtoms() const = 0;
	virtual void storeAtomLocations() = 0;

  protected:

	Minimizer(const mmcif::Polymer& poly, int first, int last,
		const mmcif::BondMap& bm, const XMap& xMap, float mapWeight,
		float plane5AtomsESD);

	double score(const AtomLocationProvider& loc);

	void addLinkRestraints(const mmcif::Monomer& a, const mmcif::Monomer& b, const std::string& linkName)
	{
		addLinkRestraints(a, b, mmcif::Link::create(linkName));
	}
	
	void addLinkRestraints(const mmcif::Monomer& a, const mmcif::Monomer& b, const mmcif::Link& link);

	template<typename R>
	double rmsz(const AtomLocationProvider& atoms, const std::vector<R>& a) const
	{
		double result = 0;
	
		if (not a.empty())
		{
			double sumZ = accumulate(a.begin(), a.end(),
				0.0, [&atoms](double sum, const R& r) { double z = r.f(atoms); return sum + z; });
			
			result = sqrt(sumZ / a.size());
		}
		
		return result;
	}

	AtomRef ref(const mmcif::Atom& atom);

	bool mElectronScattering = false;	// TODO: use!

	std::vector<mmcif::Atom> mAtoms, mReferencedAtoms;
	std::vector<size_t> mRef2AtomIndex;
	std::map<std::string,AtomRef> mRefIndex;
	
	std::vector<BondRestraint> mBondRestraints;
	std::vector<AngleRestraint> mAngleRestraints;
	std::vector<TorsionRestraint> mTorsionRestraints;
	std::vector<TransPeptideRestraint> mTransPeptideRestraints;
	std::vector<ChiralVolumeRestraint> mChiralVolumeRestraints;
	std::vector<PlanarityRestraint> mPlanarityRestraints;
	std::vector<NonBondedContactRestraint> mNonBondedContactRestraints;
	std::unique_ptr<DensityRestraint> mDensityRestraint;
	
	std::vector<Restraint*> mRestraints;
};

