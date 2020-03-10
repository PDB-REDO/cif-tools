#pragma once

#include <set>

#include <zeep/xml/serialize.hpp>

#include "cif++/Structure.h"

class NotAHBondSet
{
  public:
	~NotAHBondSet() {}

	NotAHBondSet(const NotAHBondSet&) = delete;
	NotAHBondSet& operator=(const NotAHBondSet&) = delete;

	static NotAHBondSet* Create();

	void Save(zeep::xml::element& inNode);
	static NotAHBondSet* Load(zeep::xml::element& inNode);
	
	bool IsHBondDonorOrAcceptor(const std::string& inMonomer,
		const std::string& inAtomID) const
	{
		return mData.count({inMonomer, inAtomID}) > 0;
	}
  
	bool operator()(const mmcif::Atom& atom) const
	{
		return atom.type() == mmcif::N and
			IsHBondDonorOrAcceptor(atom.comp().id(), atom.labelAtomId());
	}
  
  private:

	struct NotAHBondAtom
	{
		std::string	monomer;
		std::string	atomId;
		
		bool operator<(const NotAHBondAtom& rhs) const
		{
			int d = monomer.compare(rhs.monomer);
			if (d == 0)
				d = atomId.compare(rhs.atomId);
			return d < 0;
		}
	};
	
	typedef std::set<NotAHBondAtom> NotAHBondAtomSet;

	NotAHBondSet(NotAHBondAtomSet&& data)
		: mData(std::move(data))
	{
	}
	
	NotAHBondAtomSet mData;
};
