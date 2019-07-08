#include "pdb-redo.h"

#include <set>
#include <regex>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include "cif++/Cif++.h"

#include "HBondTraits.h"

using namespace std;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;

NotAHBondSet* NotAHBondSet::Create()
{
	const char* clibdMon = getenv("CLIBD_MON");
	if (clibdMon == nullptr)
		throw runtime_error("Cannot locate peptide list, please source the CCP4 environment");
	
	fs::path dir(clibdMon);
	
	const regex kNotAHBondEnergyTypeRX("NR(?:5|55|56|6|66)");
	NotAHBondAtomSet atomSet;
	
	for (auto d = fs::directory_iterator(dir); d != fs::directory_iterator(); ++d)
	{
		if (d->path().leaf().string().length() != 1)
			continue;
		
		for (auto e = fs::directory_iterator(d->path()); e != fs::directory_iterator(); ++e)
		{
			if (e->path().extension().string() != ".cif")
				continue;
			
			fs::ifstream file(e->path());
			if (not file.is_open())
				throw runtime_error("Could not open file: " + e->path().string());
			
			cif::File data(file);
			
			for (auto& compName: data["comp_list"]["chem_comp"])
			{
				string compId = compName["id"].as<string>();
				
				try
				{
					for (auto& compAtom: data["comp_" + compId]["chem_comp_atom"])
					{
						string atomType = compAtom["type_symbol"].as<string>();
						string typeEnergy = compAtom["type_energy"].as<string>();
					
						if (atomType == "N" and (typeEnergy == "N" or regex_match(typeEnergy, kNotAHBondEnergyTypeRX)))
							atomSet.insert({ compId, compAtom["atom_id"].as<string>() });
					}
				}
				catch (const exception& ex)
				{
					cerr << "Error reading " << e->path() << ": " << ex.what() << endl;
				}
			}
		}
	}
	
	return new NotAHBondSet(move(atomSet));
}

void NotAHBondSet::Save(zeep::xml::container& inNode)
{
	using zeep::xml::element;
	
	auto setNode = new element("not-a-hbond");
	
	for (auto& atom: mData)
	{
		auto atomNode = new element("atom");
		atomNode->set_attribute("monomer", atom.monomer);
		atomNode->set_attribute("atom_id", atom.atomId);
		setNode->append(atomNode);
	}
	
	inNode.append(setNode);
}

NotAHBondSet* NotAHBondSet::Load(zeep::xml::container& inNode)
{
	NotAHBondAtomSet data;
	
	auto setNode = inNode.find_first("not-a-hbond");
	if (setNode != nullptr)
	{
		for (auto atom: *setNode)
		{
			if (atom->name() != "atom")
				continue;
			data.insert({ atom->get_attribute("monomer"), atom->get_attribute("atom_id") });
		}
	}
	
	return new NotAHBondSet(move(data));
}
