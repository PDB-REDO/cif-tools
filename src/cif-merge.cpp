/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "cif-tools.hpp"

#include <sys/wait.h>

#include <fstream>
#include <functional>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <cif++.hpp>

namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;

using cif::iequals;

// --------------------------------------------------------------------

void updateEntryID(cif::file& target, const std::string& entryID)
{
	auto& db = target.front();
	
	if (db.name() != entryID)
		db.set_name(entryID);

	for (auto r: db["entry"])
		r["id"] = entryID;
}

void transplant(cif::file& target, cif::file& donor)
{
	auto& dbt = target.front();
	auto& dbd = donor.front();

	for (auto c: { "struct", "struct_keywords", "audit_author", "citation", "citation_author", "diffrn", "diffrn_radiation", "diffrn_radiation_wavelength" })
	{
		auto cd = dbd.get(c);
		if (cd == nullptr or cd->empty())
			continue;

		auto ct = dbt.get(c);
		if (ct == nullptr)
		{
			dbt.emplace(c);
			ct = dbt.get(c);
		}
		
		ct->clear();
		
		for (auto r: *cd)
			ct->emplace(r);
	}
	
	std::string exptlMethod;
	auto dExplt = dbd["exptl"].find(cif::key("entry_id") == dbd.name());
	if (dExplt.size() != 1)
		throw std::runtime_error("Invalid number of exptl records in donor file, should be exactly one this version of cif-merge");
	cif::tie(exptlMethod) = dExplt.front().get("method");
	auto tExplt = dbt["exptl"].find(cif::key("entry_id") == dbt.name());
	if (tExplt.empty())
	{
		auto c = dbt.emplace("exptl");
		std::get<0>(c)->emplace({
			{ "entry_id", dbt.name() },
			{ "method", exptlMethod }
		});
	}
	else
		tExplt.front()["method"] = exptlMethod;
	
	// create a mapping for the entity_ids in both files

	const std::map<std::string,const char*> kSrcMap{
		{ "man", "entity_src_gen" },
		{ "nat", "entity_src_nat" },
		{ "syn", "pdbx_entity_src_syn" }
	};
	
	std::map<std::string,std::string> d2tEntityIds;
	auto& targetEntity = dbt["entity"];
	
	for (auto r : targetEntity)
	{
		std::string id, type, dEntityID;
		cif::tie(id, type) = r.get("id", "type");
		
		if (iequals(type, "polymer"))
		{
			auto t = dbt["entity_poly"].find1(cif::key("entity_id") == id);
			
			std::string polyType, seq;
			cif::tie(polyType, seq) = t.get("type", "pdbx_seq_one_letter_code");
			
			auto d = dbd["entity_poly"].find1(cif::key("type") == polyType and cif::key("pdbx_seq_one_letter_code") == seq);
			
			if (d.empty())
			{
				if (cif::VERBOSE)
					std::cerr << "Cannot map entity " << id << " in target file to an entity in the donor" << std::endl;
				continue;
			}
			
			dEntityID = d["entity_id"].as<std::string>();
			
			// copy over refseq
			
			auto sr = dbd["struct_ref"].find1(cif::key("entity_id") == dEntityID);
			if (not sr.empty())
			{
				sr["entity_id"] = id;
				dbt["struct_ref"].emplace(sr);
				
				std::string refID = sr["id"].as<std::string>();
				
				for (auto r: dbd["struct_ref_seq"].find(cif::key("ref_id") == refID))
					dbt["struct_ref_seq"].emplace(r);
			}
		}
		else if (iequals(type, "non-polymer"))
		{
			auto t = dbt["pdbx_entity_nonpoly"].find1(cif::key("entity_id") == id);
			
			std::string compID;
			cif::tie(compID) = t.get("comp_id");
			
			auto d = dbd["pdbx_entity_nonpoly"].find1(cif::key("comp_id") == compID);
			
			if (d.empty())
			{
				if (cif::VERBOSE)
					std::cerr << "Cannot map entity " << id << " in target file to an entity in the donor" << std::endl;
				continue;
			}
			
			cif::tie(dEntityID) = d.get("entity_id");
		}
		else if (iequals(type, "water"))
		{
			dEntityID = dbd["entity"].find1<std::string>(cif::key("type") == type, "id");
		}
		else if (cif::VERBOSE)
			std::cerr << "Unsupported entity type: " << type << std::endl;
		
		if (dEntityID.empty())
			continue;

		const auto &[srcMethod, description, weight] =
			dbd["entity"].find1<std::string,std::string,std::string>(cif::key("id") == dEntityID, "src_method", "pdbx_description", "formula_weight");
		
		r["src_method"] = srcMethod;
		r["pdbx_description"] = description;
		r["formula_weight"] = weight;
		
		if (kSrcMap.count(srcMethod))
		{
			std::string srcRec = kSrcMap.at(srcMethod);

			auto d = dbd[srcRec].find1(cif::key("entity_id") == dEntityID);
			if (not d.empty())
			{
				d["entity_id"] = id;
				dbt.emplace(srcRec);
				dbt[srcRec].emplace(d);
			}
		}
	}
}

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-merge [options] inputFile donorFile ");
	visible_options.add_options()
		("help,h",									"Display help message")
		("version",									"Print version")
		("verbose,v",								"Verbose output")
		("input,i",		po::value<std::string>(),	"Modified PDB file")
		("output,o",	po::value<std::string>(),	"Output file, default is stdout (terminal)")
		("donor",		po::value<std::string>(),	"CIF file (or PDB ID for this file) containing the data to collect data from")
		// ("dict",		po::value<std::string>(),	"Dictionary file containing restraints for residues in this specific target")
		;
	

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("donor", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		write_version_string(std::cout, vm.count("verbose"));
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm.count("donor") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// Load dict, if any
	
	// if (vm.count("dict"))
	// 	c::CompoundFactory::instance().pushDictionary(vm["dict"].as<std::string>());

	// Read input file
	cif::file cf{vm["input"].as<std::string>()};
	
	// Read donor file
	cif::file df{vm["donor"].as<std::string>()};

	updateEntryID(cf, df.front().name());
	transplant(cf, df);
	
	if (vm.count("output"))
		cf.save(vm["output"].as<std::string>());
	else
		cf.save(std::cout);

	return 0;	
}
