/* 
	Created by: Maarten L. Hekkelman
	Date: maandag 07 januari, 2019

	Skip list for e.g. pepflip
*/

#pragma once

#include <boost/filesystem/path.hpp>

#include "cif++/Structure.h"

struct SkipResidue
{
	std::string	asymID;
	int			seqID;
	std::string	monID;
};

typedef std::vector<SkipResidue> SkipList;

enum SkipListNumberingScheme
{
	sl_Label, sl_PDB //, sl_Auth
};

// --------------------------------------------------------------------

SkipList readSkipList(const boost::filesystem::path& file,
	const mmcif::Structure& structure);

void writeSkipList(const boost::filesystem::path& file,
	const mmcif::Structure& structure, const SkipList& list,
	SkipListNumberingScheme scheme = sl_Label);
