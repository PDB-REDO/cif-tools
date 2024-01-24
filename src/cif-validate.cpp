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

#include <fstream>

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include "revision.hpp"

class dummy_parser : public cif::sac_parser
{
  public:
	dummy_parser(std::istream &is)
		: sac_parser(is)
	{
	}

	void produce_datablock(std::string_view name) override
	{
	}

	void produce_category(std::string_view name) override
	{
	}

	void produce_row() override
	{
	}

	void produce_item(std::string_view category, std::string_view item, std::string_view value) override
	{
	}
};

void validate_syntax(const std::string &file)
{
	cif::gzio::ifstream in(file);
	dummy_parser parser(in);
	parser.parse_file();
}

int pr_main(int argc, char *argv[])
{
	auto &config = mcfp::config::instance();

	config.init("cif-validate [options...] file [output-file]",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option<std::string>("dict", "mmcif_pdbx.dic", "The mmCIF dictionary to use, can be either mmcif_ddl, mmcif_pdbx or a path to the actual dictionary file"),
		mcfp::make_option("validate-pdbx", "Validate PDBx categories, assuming mmcif_pdbx, and checks entity, entity_poly, entity_poly_seq and pdbx_poly_seq_scheme for consistency with atom_site"),
		mcfp::make_option("validate-links", "Validate all links"),
		mcfp::make_option("validate-data-comp", "Default is to skip data_comp_XXX datablocks, use this flag to force validation"),
		mcfp::make_option("syntax-only", "Quickly check to see if the syntax is correct"),
		mcfp::make_option("verbose,v", "Verbose output, repeat to increase verbosity level"),
		mcfp::make_option("print", "Print the reformatted file, to stdout or, when specified, to 'output-file'"));

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().empty())
	{
		std::cerr << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	cif::VERBOSE = config.count("verbose");

	int result = 0;

	if (config.has("syntax-only") and not config.has("print"))
		validate_syntax(config.operands().front());
	else
	{
		cif::file f(config.operands().front());

		if (not config.has("syntax-only"))
		{
			if (config.count("dict"))
				f.load_dictionary(config.get<std::string>("dict"));
			else if (config.has("validate-pdbx"))
				f.load_dictionary("mmcif_pdbx");

			if (f.get_validator() == nullptr)
			{
				std::cerr << "No validator, please specify a dictionary to use using the --dict option" << std::endl
						  << "However, the syntax seems to be OK" << std::endl;
			}
			else
			{
				for (auto &db : f)
				{
					if (cif::starts_with(db.name(), "comp_") and not config.has("validate-data-comp"))
						continue;

					if (not db.is_valid())
						result = 1;;
				}

				if (config.has("validate-links"))
					f.validate_links();
				
				if (config.has("validate-pdbx"))
					result = result and cif::pdb::is_valid_pdbx_file(f);
			}
		}

		if (config.has("print"))
		{
			if (config.operands().size() == 1)
				f.save(std::cout);
			else
			{
				std::ofstream out(config.operands()[1]);
				if (not out.is_open())
					std::cerr << "Could not open output file" << std::endl;
				else
					f.save(out);
			}
		}
	}

	return result;
}
