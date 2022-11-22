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

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

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
	if (target.empty() or donor.empty())
		throw std::runtime_error("empty files?");

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
				if (cif::VERBOSE > 0)
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
				if (cif::VERBOSE > 0)
					std::cerr << "Cannot map entity " << id << " in target file to an entity in the donor" << std::endl;
				continue;
			}
			
			cif::tie(dEntityID) = d.get("entity_id");
		}
		else if (iequals(type, "water"))
		{
			dEntityID = dbd["entity"].find1<std::string>(cif::key("type") == type, "id");
		}
		else if (cif::VERBOSE > 0)
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
	auto &config = mcfp::config::instance();

	config.init(
		"cif-merge input-file donor-file [output-file]",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "Verbose output"),
		mcfp::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)")
	);

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().size() < 2)
	{
		std::cerr << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");
	if (config.has("debug"))
		cif::VERBOSE = config.get<int>("debug");

	// Load dict, if any
	
	// if (vm.count("dict"))
	// 	c::CompoundFactory::instance().pushDictionary(vm["dict"].as<std::string>());

	// Read input file
	cif::gzio::ifstream in(config.operands().front());
	if (not in.is_open())
		throw std::runtime_error("Could not open input file");
	
	cif::gzio::ifstream donor(config.operands()[1]);
	if (not donor.is_open())
		throw std::runtime_error("Could not open donor file");
	
	cif::file cf{in};
	cif::file df{donor};

	updateEntryID(cf, df.front().name());
	transplant(cf, df);
	
	if (config.operands().size() == 3)
		cf.save(config.operands().back());
	else
		cf.save(std::cout);

	return 0;	
}
