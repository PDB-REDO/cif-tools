/* include/cif++/Config.hpp.  Generated from Config.hpp.in by configure.  */
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

#include <sys/time.h>
#include <sys/resource.h>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <regex>

#include "cif++/Cif++.hpp"
#include "cif++/CifUtils.hpp"

string VERSION_STRING;

int pr_main(int argc, char* argv[]);

// --------------------------------------------------------------------

ostream& operator<<(ostream& os, const struct timeval& t)
{
	uint64_t s = t.tv_sec;
	if (s > 24 * 60 * 60)
	{
		uint32_t days = s / (24 * 60 * 60);
		os << days << "d ";
		s %= 24 * 60 * 60;
	}
	
	if (s > 60 * 60)
	{
		uint32_t hours = s / (60 * 60);
		os << hours << "h ";
		s %= 60 * 60;
	}
	
	if (s > 60)
	{
		uint32_t minutes = s / 60;
		os << minutes << "m ";
		s %= 60;
	}
	
	double ss = s + 1e-6 * t.tv_usec;
	
	os << fixed << setprecision(1) << ss << 's';

	return os;
}

ostream& operator<<(ostream& os, const chrono::duration<double>& t)
{
	uint64_t s = static_cast<uint64_t>(trunc(t.count()));
	if (s > 24 * 60 * 60)
	{
		uint32_t days = s / (24 * 60 * 60);
		os << days << "d ";
		s %= 24 * 60 * 60;
	}
	
	if (s > 60 * 60)
	{
		uint32_t hours = s / (60 * 60);
		os << hours << "h ";
		s %= 60 * 60;
	}
	
	if (s > 60)
	{
		uint32_t minutes = s / 60;
		os << minutes << "m ";
		s %= 60;
	}
	
	double ss = s + 1e-6 * (t.count() - s);
	
	os << fixed << setprecision(1) << ss << 's';

	return os;
}

class RUsage
{
  public:
	~RUsage()
	{
		if (cif::VERBOSE)
		{
			struct rusage u;
			auto end = std::chrono::system_clock::now();
			chrono::duration<double> diff = end - start;
			
			if (getrusage(RUSAGE_SELF, &u) == 0)
				cerr << "CPU usage: "
					<< u.ru_utime << " user, "
					<< u.ru_stime << " system, "
					<< diff << " wall" << endl;
			else
				perror("Failed to get rusage");
		}
	}

	chrono::time_point<chrono::system_clock>	start = std::chrono::system_clock::now();
};

// --------------------------------------------------------------------

namespace {
	std::string gVersionNr, gVersionDate;
}

void load_version_info()
{
	const regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)"),
		rxVersionDate(R"(Date: +(\d{4}-\d{2}-\d{2}).*)");

#if USE_RSRC

	auto version = cif::rsrc_loader::load("version.txt");
	if (not version)
	{
		VERSION_STRING = "unknown version, version resource is missing";
		return;
	}

	struct membuf : public streambuf
	{
		membuf(char* data, size_t length)       { this->setg(data, data, data + length); }
	} buffer(const_cast<char*>(version.data()), version.size());

#else

#include "revision.hpp"

	struct membuf : public streambuf
	{
		membuf(char* data, size_t length)       { this->setg(data, data, data + length); }
	} buffer(const_cast<char*>(kRevision), sizeof(kRevision));

#endif
	istream is(&buffer);

	string line;

	while (getline(is, line))
	{
		smatch m;

		if (regex_match(line, m, rxVersionNr))
		{
			gVersionNr = m[1];
			if (m[2].matched)
				gVersionNr += '*';
			continue;
		}

		if (regex_match(line, m, rxVersionDate))
		{
			gVersionDate = m[1];
			continue;
		}
	}

	if (not VERSION_STRING.empty())
		VERSION_STRING += "\n";
	VERSION_STRING += gVersionNr + " " + gVersionDate;
}

string get_version_nr()
{
	return gVersionNr;
}

string get_version_date()
{
	return gVersionDate;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what (const exception& e)
{
	cerr << e.what() << endl;
	try
	{
		rethrow_if_nested(e);
	}
	catch (const exception& nested)
	{
		cerr << " >> ";
		print_what(nested);
	}
}

int main(int argc, char* argv[])
{
	int result = -1;
	
	RUsage r;
	
	try
	{
		cif::rsrc_loader::init({
			{ cif::rsrc_loader_type::file, "." },
#if defined DATADIR
			{ cif::rsrc_loader_type::file, DATADIR },
#endif
#if USE_RSRC
			{ cif::rsrc_loader_type::mrsrc, "", { gResourceIndex, gResourceData, gResourceName } }
#endif
		});

		load_version_info();
		
		result = pr_main(argc, argv);
	}
	catch (exception& ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
