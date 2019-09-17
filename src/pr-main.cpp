#include "pdb-redo.h"

#include <sys/time.h>
#include <sys/resource.h>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <regex>

#include "cif++/mrsrc.h"

using namespace std;

string VERSION;

int pr_main(int argc, char* argv[]);

int VERBOSE;

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
		if (VERBOSE)
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

static string gVersionNr, gVersionDate;

void load_version_info()
{
	mrsrc::rsrc version("version.txt");
	if (not version)
		VERSION = "unknown version, version resource is missing";
	else
	{
		struct membuf : public streambuf
		{
			membuf(char* data, size_t length)		{ this->setg(data, data, data + length); }
		} buffer(const_cast<char*>(version.data()), version.size());
		
		istream is(&buffer);
		string line;
		regex
			rxVersionNr(R"(Last Changed Rev: (\d+))"),
			rxVersionDate(R"(Last Changed Date: (\d{4}-\d{2}-\d{2}).*)");

		while (getline(is, line))
		{
			smatch m;

			if (regex_match(line, m, rxVersionNr))
			{
				gVersionNr = m[1];
				continue;
			}

			if (regex_match(line, m, rxVersionDate))
			{
				gVersionDate = m[1];
				continue;
			}
		}

		VERSION = gVersionNr + " " + gVersionDate;
	}
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
