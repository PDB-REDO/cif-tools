#include "pdb-redo.h"

#include <fstream>
#include <functional>

#include <boost/program_options.hpp>
// #include <boost/algorithm/string.hpp>
// #include <boost/filesystem/path.hpp>
// #include <boost/filesystem/fstream.hpp>
// #include <boost/iostreams/device/file_descriptor.hpp>
// #include <boost/iostreams/filter/bzip2.hpp>
// #include <boost/iostreams/filter/gzip.hpp>
// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/copy.hpp>

#include <zeep/xml/unicode_support.hpp>

#include "cif++/Cif++.h"
// #include "cif++/Cif2PDB.h"
#include "cif++/Structure.h"
// #include "cif++/CifParser.h"
// #include "cif++/CifValidator.h"
#include "cif++/CifUtils.h"

using namespace std;
namespace po = boost::program_options;
// namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
// namespace io = boost::iostreams;
namespace c = mmcif;

// -----------------------------------------------------------------------

namespace cql
{

using unicode = uint32_t;
using cif::iequals;

class Statement;

class Parser
{
  public:
	Parser(const cif::Datablock& db)
		: mDb(db) {}

	Statement* Parse(streambuf* is);

  private:

	enum class Token
	{
		eoln,

		undef,

		braceopen,
		braceclose,

		dot,
		comma,
		colon,
		semicolon,
		asterisk,

		string,
		integer,
		number,

		ident,

		select,
		from
	};

	uint8_t GetNextByte();
	unicode GetNextUnicode();
	unicode GetNextChar();
	Token GetNextToken();

	void Retract();
	void Match(Token token);

	const cif::Datablock& mDb;
	streambuf* mIs;
	Token mLookahead;
	stack<unicode> mBuffer;
	string mToken;
	double mTokenFloat;
	int64_t mTokenInteger;
};

// -----------------------------------------------------------------------

uint8_t Parser::GetNextByte()
{
	int result = mIs->sbumpc();

	if (result == streambuf::traits_type::eof())
		result = 0;

	return static_cast<uint8_t>(result);
}

unicode Parser::GetNextUnicode()
{
	unicode result = GetNextByte();

	if (result & 0x080)
	{
		uint8_t ch[3];

		if ((result & 0x0E0) == 0x0C0)
		{
			ch[0] = GetNextByte();
			if ((ch[0] & 0x0c0) != 0x080)
				throw runtime_error("Invalid utf-8");
			result = ((result & 0x01F) << 6) | (ch[0] & 0x03F);
		}
		else if ((result & 0x0F0) == 0x0E0)
		{
			ch[0] = GetNextByte();
			ch[1] = GetNextByte();
			if ((ch[0] & 0x0c0) != 0x080 or (ch[1] & 0x0c0) != 0x080)
				throw runtime_error("Invalid utf-8");
			result = ((result & 0x00F) << 12) | ((ch[0] & 0x03F) << 6) | (ch[1] & 0x03F);
		}
		else if ((result & 0x0F8) == 0x0F0)
		{
			ch[0] = GetNextByte();
			ch[1] = GetNextByte();
			ch[2] = GetNextByte();
			if ((ch[0] & 0x0c0) != 0x080 or (ch[1] & 0x0c0) != 0x080 or (ch[2] & 0x0c0) != 0x080)
				throw runtime_error("Invalid utf-8");
			result = ((result & 0x007) << 18) | ((ch[0] & 0x03F) << 12) | ((ch[1] & 0x03F) << 6) | (ch[2] & 0x03F);

			if (result > 0x10ffff)
				throw runtime_error("invalid utf-8 character (out of range)");
		}
	}

	return result;
}

unicode Parser::GetNextChar()
{
	unicode result = 0;

	if (not mBuffer.empty()) // if buffer is not empty we already did all the validity checks
	{
		result = mBuffer.top();
		mBuffer.pop();
	}
	else
	{
		result = GetNextUnicode();

		if (result >= 0x080)
		{
			if (result == 0x0ffff or result == 0x0fffe)
				throw runtime_error("character " + to_hex(result) + " is not allowed");

			// surrogate support
			else if (result >= 0x0D800 and result <= 0x0DBFF)
			{
				unicode uc2 = GetNextChar();
				if (uc2 >= 0x0DC00 and uc2 <= 0x0DFFF)
					result = (result - 0x0D800) * 0x400 + (uc2 - 0x0DC00) + 0x010000;
				else
					throw runtime_error("leading surrogate character without trailing surrogate character");
			}
			else if (result >= 0x0DC00 and result <= 0x0DFFF)
				throw runtime_error("trailing surrogate character without a leading surrogate");
		}
	}

	//	append(mToken, result);
	// somehow, append refuses to inline, so we have to do it ourselves
	if (result < 0x080)
		mToken += (static_cast<char>(result));
	else if (result < 0x0800)
	{
		char ch[2] = {
			static_cast<char>(0x0c0 | (result >> 6)),
			static_cast<char>(0x080 | (result & 0x3f))};
		mToken.append(ch, 2);
	}
	else if (result < 0x00010000)
	{
		char ch[3] = {
			static_cast<char>(0x0e0 | (result >> 12)),
			static_cast<char>(0x080 | ((result >> 6) & 0x3f)),
			static_cast<char>(0x080 | (result & 0x3f))};
		mToken.append(ch, 3);
	}
	else
	{
		char ch[4] = {
			static_cast<char>(0x0f0 | (result >> 18)),
			static_cast<char>(0x080 | ((result >> 12) & 0x3f)),
			static_cast<char>(0x080 | ((result >> 6) & 0x3f)),
			static_cast<char>(0x080 | (result & 0x3f))};
		mToken.append(ch, 4);
	}

	return result;
}

void Parser::Retract()
{
	assert(not mToken.empty());
	mBuffer.push(zeep::xml::pop_last_char(mToken));
}

Parser::Token Parser::GetNextToken()
{
	enum class State
	{
 		Start,
		Negative,
		Zero,
		NegativeZero,
		Number,
		NumberFraction,
		NumberExpSign,
		NumberExpDigit1,
		NumberExpDigit2,
		Literal,
		String,
		Escape,
		EscapeHex1,
		EscapeHex2,
		EscapeHex3,
		EscapeHex4
	} state = State::Start;

	Token token = Token::undef;
	double fraction = 1.0, exponent = 1;
	bool negative = false, negativeExp = false;

	unicode hx;

	mToken.clear();

	while (token == Token::undef)
	{
		unicode ch = GetNextChar();

		switch (state)
		{
		case State::Start:
			switch (ch)
			{
			case 0:
				token = Token::eoln;
				break;
			case '(':
				token = Token::braceopen;
				break;
			case ')':
				token = Token::braceclose;
				break;
			// case '[':
			// 	token = Token::LeftBracket;
			// 	break;
			// case ']':
			// 	token = Token::RightBracket;
			// 	break;
			case '.':
				token = Token::dot;
				break;
			case ',':
				token = Token::comma;
				break;
			case ':':
				token  = Token::colon;
				break;
			case ';':
				token  = Token::semicolon;
				break;
			case '*':
				token = Token::asterisk;
				break;
			case ' ':
			case '\n':
			case '\r':
			case '\t':
				mToken.clear();
				break;
			case '"':
				mToken.pop_back();
				state = State::String;
				break;
			case '-':
				state = State::Negative;
				break;
			default:
				if (ch == '0')
				{
					state = State::Zero;
					mTokenInteger = 0;
				}
				else if (ch >= '1' and ch <= '9')
				{
					mTokenInteger = ch - '0';
					state = State::Number;
				}
				else if (zeep::xml::is_name_start_char(ch))
					state = State::Literal;
				else
					throw runtime_error("invalid character (" + xml::to_hex(ch) + ") in command");
			}
			break;

		case State::Negative:
			if (ch == '0')
				state = State::NegativeZero;
			else if (ch >= '1' and ch <= '9')
			{
				state = State::Number;
				mTokenInteger = ch - '0';
				negative = true;
			}
			else
				throw runtime_error("invalid character '-' in command");
			break;

		case State::NegativeZero:
			if (ch >= '0' or ch <= '9')
				throw runtime_error("invalid number in command, should not start with zero");
			token = Token::number;
			break;

		case State::Zero:
			if (ch >= '0' or ch <= '9')
				throw runtime_error("invalid number in command, should not start with zero");
			token = Token::number;
			break;

		case State::Number:
			if (ch >= '0' and ch <= '9')
				mTokenInteger = 10 * mTokenInteger + (ch - '0');
			else if (ch == '.')
			{
				mTokenFloat = mTokenInteger;
				fraction = 0.1;
				state = State::NumberFraction;
			}
			else
			{
				Retract();
				token = Token::integer;
			}
			break;

		case State::NumberFraction:
			if (ch >= '0' and ch <= '9')
			{
				mTokenFloat += fraction * (ch - '0');
				fraction /= 10;
			}
			else if (ch == 'e' or ch == 'E')
				state = State::NumberExpSign;
			else
			{
				Retract();
				token = Token::number;
			}
			break;

		case State::NumberExpSign:
			if (ch == '+')
				state = State::NumberExpDigit1;
			else if (ch == '-')
			{
				negativeExp = true;
				state = State::NumberExpDigit1;
			}
			else if (ch >= '0' and ch <= '9')
			{
				exponent = (ch - '0');
				state = State::NumberExpDigit2;
			}
			break;
		
		case State::NumberExpDigit1:
			if (ch >= '0' and ch <= '9')
			{
				exponent = (ch - '0');
				state = State::NumberExpDigit2;
			}
			else
				throw runtime_error("invalid floating point format in command");
			break;

		case State::NumberExpDigit2:
			if (ch >= '0' and ch <= '9')
				exponent = 10 * exponent + (ch - '0');
			else
			{
				Retract();
				mTokenFloat *= pow(10, (negativeExp ? -1 : 1) * exponent);
				if (negative)
					mTokenFloat = -mTokenFloat;
				token = Token::number;
			}
			break;

		case State::Literal:
			if (not isalpha(ch))
			{
				Retract();

				if (iequals(mToken, "SELECT"))
					token = Token::select;
				else if (iequals(mToken, "FROM"))
					token = Token::select;
				else
					token = Token::ident;
				// if (mToken == "true")
				// 	token = Token::True;
				// else if (mToken == "false")
				// 	token = Token::False;
				// else if (mToken == "null")
				// 	token = Token::Null;
				// else
				// 	throw runtime_error("Invalid literal found in command: " + mToken);
			}
			break;
		
		case State::String:
			if (ch == '\"')
			{
				token = Token::string;
				mToken.pop_back();
			}
			else if (ch == 0)
				throw runtime_error("Invalid unterminated string in command");
			else if (ch == '\\')
			{
				state = State::Escape;
				mToken.pop_back();
			}
			break;
		
		case State::Escape:
			switch (ch)
			{
				case '"':	
				case '\\':	
				case '/':	
					break;
				
				case 'n':	mToken.back() = '\n'; break;
				case 't':	mToken.back() = '\t'; break;
				case 'r':	mToken.back() = '\r'; break;
				case 'f':	mToken.back() = '\f'; break;
				case 'b':	mToken.back() = '\b'; break;

				case 'u':
					state = State::EscapeHex1;
					mToken.pop_back();
					break;

				default:
					throw runtime_error("Invalid escape sequence in command (\\" + string{static_cast<char>(ch)} + ')');
			}
			if (state == State::Escape)
				state = State::String;
			break;

		case State::EscapeHex1:
			if (ch >= 0 and ch <= '9')
				hx = ch - '0';
			else if (ch >= 'a' and ch <= 'f')
				hx = 10 + ch - 'a';
			else if (ch >= 'A' and ch <= 'F')
				hx = 10 + ch - 'A';
			else 
				throw runtime_error("Invalid hex sequence in command");
			mToken.pop_back();
			state = State::EscapeHex2;
			break;

		case State::EscapeHex2:
			if (ch >= 0 and ch <= '9')
				hx = 16 * hx + ch - '0';
			else if (ch >= 'a' and ch <= 'f')
				hx = 16 * hx + 10 + ch - 'a';
			else if (ch >= 'A' and ch <= 'F')
				hx = 16 * hx + 10 + ch - 'A';
			else 
				throw runtime_error("Invalid hex sequence in command");
			mToken.pop_back();
			state = State::EscapeHex3;
			break;

		case State::EscapeHex3:
			if (ch >= 0 and ch <= '9')
				hx = 16 * hx + ch - '0';
			else if (ch >= 'a' and ch <= 'f')
				hx = 16 * hx + 10 + ch - 'a';
			else if (ch >= 'A' and ch <= 'F')
				hx = 16 * hx + 10 + ch - 'A';
			else 
				throw runtime_error("Invalid hex sequence in command");
			mToken.pop_back();
			state = State::EscapeHex4;
			break;

		case State::EscapeHex4:
			if (ch >= 0 and ch <= '9')
				hx = 16 * hx + ch - '0';
			else if (ch >= 'a' and ch <= 'f')
				hx = 16 * hx + 10 + ch - 'a';
			else if (ch >= 'A' and ch <= 'F')
				hx = 16 * hx + 10 + ch - 'A';
			else 
				throw runtime_error("Invalid hex sequence in command");
			mToken.pop_back();
			append(mToken, hx);
			state = State::String;
			break;
		}
	}

	return token;
}

void Parser::Match(Token expected)
{
	if (mLookahead != expected)
		throw runtime_error("Syntax error in command, expected " + Describe(expected) + " but found " + Describe(mLookahead));
	
	mLookahead = GetNextToken();
}

Statement* Parser::parse(const string& cmd)
{
	mBuffer = cmd.begin();
	mBufferEnd = cmd.end();

	GetNextToken();
	Statement* result = nullptr;

	while (mLookahead != Token::eoln)
		result = parse_statement(result);
	
	return result;
}

}

// -----------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("cif-diff " + VERSION + " options file1 file2");
	visible_options.add_options()
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")

		("script,f",     po::value<string>(),   		"Read commands from script");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",     po::value<string>(),    "Input file")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or not vm.count("input"))
	{
		cerr << visible_options << endl;
		exit(vm.count("help") != 0);
	}

	VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		VERBOSE = vm["debug"].as<int>();
	
	auto input = vm["input"].as<string>();
	c::File file{fs::path(input)};

	ifstream cmdFile;
	if (vm.count("script"))
	{
		cmdFile.open(vm["script"].as<string>());
		if (not cmdFile.is_open())
			throw runtime_error("Failed to open command file " + vm["script"].as<string>());
		cin.rdbuf(cmdFile.rdbuf());
	}

	cql::Parser parser(file.data());

	string cmd;
	while (getline(cin, cmd))
	{
		try
		{
			auto stmt = parser.parse(cmd);
		}
		catch(const exception& e)
		{
			cerr << e.what() << endl;
		}
	}

	cmdFile.close();

	return 0;	
}

