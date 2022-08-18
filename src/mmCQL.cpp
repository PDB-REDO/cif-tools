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

#include <fstream>
#include <filesystem>
#include <functional>
#include <unordered_set>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "cif++.hpp"

namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;

using unicode = char32_t;

namespace zeep {

// inlines

/// \brief Append a single unicode character to an utf-8 string
inline void append(std::string& s, unicode uc)
{
	if (uc < 0x080)
		s += (static_cast<char>(uc));
	else if (uc < 0x0800)
	{
		char ch[2] = {
			static_cast<char>(0x0c0 | (uc >> 6)),
			static_cast<char>(0x080 | (uc & 0x3f))
		};
		s.append(ch, 2);
	}
	else if (uc < 0x00010000)
	{
		char ch[3] = {
			static_cast<char>(0x0e0 | (uc >> 12)),
			static_cast<char>(0x080 | ((uc >> 6) & 0x3f)),
			static_cast<char>(0x080 | (uc & 0x3f))
		};
		s.append(ch, 3);
	}
	else
	{
		char ch[4] = {
			static_cast<char>(0x0f0 | (uc >> 18)),
			static_cast<char>(0x080 | ((uc >> 12) & 0x3f)),
			static_cast<char>(0x080 | ((uc >> 6) & 0x3f)),
			static_cast<char>(0x080 | (uc & 0x3f))
		};
		s.append(ch, 4);
	}
}

/// \brief remove the last unicode character from an utf-8 string
inline unicode pop_last_char(std::string& s)
{
	unicode result = 0;

	if (not s.empty())
	{
		std::string::iterator ch = s.end() - 1;
		
		if ((*ch & 0x0080) == 0)
		{
			result = *ch;
			s.erase(ch);
		}
		else
		{
			int o = 0;
			
			do
			{
				result |= (*ch & 0x03F) << o;
				o += 6;
				--ch;
			}
			while (ch != s.begin() and (*ch & 0x0C0) == 0x080);
			
			switch (o)
			{
				case  6: result |= (*ch & 0x01F) <<  6; break;
				case 12: result |= (*ch & 0x00F) << 12; break;
				case 18: result |= (*ch & 0x007) << 18; break;
			}
			
			s.erase(ch, s.end());
		}
	}
	
	return result;
}

// this code only works if the input is valid utf-8
/// \brief return the first unicode and the advanced pointer from a string
template<typename Iter>
std::tuple<unicode,Iter> get_first_char(Iter ptr)
{
	unicode result = static_cast<unsigned char>(*ptr);
	++ptr;

	if (result > 0x07f)
	{
		unsigned char ch[3];
		
		if ((result & 0x0E0) == 0x0C0)
		{
			ch[0] = static_cast<unsigned char>(*ptr); ++ptr;
			result = ((result & 0x01F) << 6) | (ch[0] & 0x03F);
		}
		else if ((result & 0x0F0) == 0x0E0)
		{
			ch[0] = static_cast<unsigned char>(*ptr); ++ptr;
			ch[1] = static_cast<unsigned char>(*ptr); ++ptr;
			result = ((result & 0x00F) << 12) | ((ch[0] & 0x03F) << 6) | (ch[1] & 0x03F);
		}
		else if ((result & 0x0F8) == 0x0F0)
		{
			ch[0] = static_cast<unsigned char>(*ptr); ++ptr;
			ch[1] = static_cast<unsigned char>(*ptr); ++ptr;
			ch[2] = static_cast<unsigned char>(*ptr); ++ptr;
			result = ((result & 0x007) << 18) | ((ch[0] & 0x03F) << 12) | ((ch[1] & 0x03F) << 6) | (ch[2] & 0x03F);
		}
	}

	return std::make_tuple(result, ptr);
}

// --------------------------------------------------------------------

inline std::string to_hex(uint32_t i)
{
	char s[sizeof(i) * 2 + 3];
	char* p = s + sizeof(s);
	*--p = 0;

	const char kHexChars[] = "0123456789abcdef";

	while (i)
	{
		*--p = kHexChars[i & 0x0F];
		i >>= 4;
	}

	*--p = 'x';
	*--p = '0';

	return p;
}	
}

// -----------------------------------------------------------------------

namespace cql
{

using unicode = uint32_t;
using cif::iequals;

class Statement;
typedef std::shared_ptr<Statement> StatementPtr;

// -----------------------------------------------------------------------

class Statement
{
  public:
	Statement(const Statement& ) = delete;
	Statement& operator=(const Statement&) = delete;

	virtual ~Statement() {}

	virtual void Execute() = 0;

  protected:
	Statement() {}
};

// -----------------------------------------------------------------------

class StatementList : public Statement
{
  public:
	StatementList() {}

	void Add(StatementPtr stmt)
	{
		mStatements.emplace_back(stmt);
	}

	virtual void Execute()
	{
		for (auto stmt: mStatements)
			stmt->Execute();
	}

  private:
	std::vector<StatementPtr> mStatements;
};

// -----------------------------------------------------------------------

class SelectStatement : public Statement
{
  public:
	SelectStatement(cif::category& category, bool distinct, std::vector<std::string>&& items, cif::condition&& where)
		: mCategory(category), mDistinct(distinct), mItems(std::move(items)), mWhere(std::move(where)) {}

	virtual void Execute()
	{
		std::vector<std::string> fields(mItems.size());
		std::unordered_set<std::string> seen;

		std::cout << ba::join(mItems, "\t") << std::endl;

		for (auto r: mCategory.find(std::move(mWhere)))
		{
			transform(mItems.begin(), mItems.end(), fields.begin(),
				[r](auto item) {
					return r[item].template as<std::string>();
					});

			std::string line = ba::join(fields, "\t");
			bool seenLine = seen.count(line);

			if (not mDistinct or not seenLine)
				std::cout << line << std::endl;

			if (mDistinct and not seenLine)
				seen.insert(line);
		}
	}

  private:
	cif::category& mCategory;
	bool mDistinct;
	std::vector<std::string> mItems;
	cif::condition mWhere;
};

// -----------------------------------------------------------------------

class DeleteStatement : public Statement
{
  public:
	DeleteStatement(cif::category& category, cif::condition&& where)
		: mCategory(category), mWhere(std::move(where)) {}

	virtual void Execute()
	{
		std::vector<cif::row_handle> remove;
		
		mWhere.prepare(mCategory);

		for (auto r: mCategory)
		{
			if (mWhere(r))
				remove.insert(remove.end(), r);
		}

		for (auto r: remove)
			mCategory.erase(r);

		std::cout << "Number of removed rows " << remove.size() << std::endl;
	}

  private:
	cif::category& mCategory;
	cif::condition mWhere;
};

// -----------------------------------------------------------------------

class UpdateStatement : public Statement
{
  public:
	UpdateStatement(cif::category& category, std::vector<std::pair<std::string,std::string>>&& itemValuePairs, cif::condition&& where)
		: mCategory(category), mItemValuePairs(std::move(itemValuePairs)), mWhere(std::move(where)) {}

	virtual void Execute()
	{
		size_t updated = 0;

		mWhere.prepare(mCategory);

		for (auto r: mCategory)
		{
			if (mWhere(r))
			{
				for (auto iv: mItemValuePairs)
					r[iv.first] = iv.second;

				++updated;
			}
		}

		std::cout << "Number of updated rows: " << updated << std::endl;
	}

  private:
	cif::category& mCategory;
	std::vector<std::pair<std::string,std::string>> mItemValuePairs;
	cif::condition mWhere;
};

// -----------------------------------------------------------------------

class Parser
{
  public:
	Parser(cif::datablock& db)
		: mDb(db) {}

	StatementPtr Parse(std::streambuf* is);

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

		eq_,
		lt_,
		le_,
		gt_,
		ge_,
		ne_,

		string,
		integer,
		number,

		ident,

		select,
		distinct,
		from,
		update,
		set,
		where,
		and_,
		or_,
		not_,
		insert,
		delete_,
		into,
		values,
		is_,
		null_
	};

	std::string Describe(Token token)
	{
		switch (token)
		{
			case Token::eoln:   	return "<EOLN>";
			case Token::undef:		return "<UNDEFINED>";
			case Token::braceopen:	return "'('";
			case Token::braceclose:	return "')'";

			case Token::dot:		return "'.'";
			case Token::comma:		return "','";
			case Token::colon:		return "':'";
			case Token::semicolon:	return "';'";
			case Token::asterisk:	return "'*'";

			case Token::eq_:		return "'='";
			case Token::lt_:		return "'<'";
			case Token::le_:		return "'<='";
			case Token::gt_:		return "'>'";
			case Token::ge_:		return "'>='";
			case Token::ne_:		return "'<>'";

			case Token::string:		return "string";
			case Token::integer:	return "integer";
			case Token::number:		return "number";
			case Token::ident:		return "identifier";

			case Token::select:		return "SELECT";
			case Token::distinct:	return "DISTINCT";
			case Token::from:		return "FROM";
			case Token::update:		return "UPDATE";
			case Token::set:		return "SET";
			case Token::where:		return "WHERE";
			case Token::and_:		return "AND";
			case Token::or_:		return "OR";
			case Token::not_:		return "NOT";
			case Token::insert:		return "INSERT";
			case Token::delete_:	return "DELETE";
			case Token::into:		return "INTO";
			case Token::values:		return "VALUES";
			case Token::is_:		return "IS";
			case Token::null_:		return "NULL";

			default:			assert(false); return "unknown token";
		}
	}

	uint8_t GetNextByte();
	unicode GetNextUnicode();
	unicode GetNextChar();
	Token GetNextToken();

	void Retract();
	void Match(Token token);

	// parser rules
	StatementPtr ParseStatement();
	StatementPtr ParseSelect();
	StatementPtr ParseDelete();
	StatementPtr ParseUpdate();
	std::vector<std::string> ParseItemList();

	cif::condition ParseWhereClause(cif::category& cat);
	cif::condition ParseNotWhereClause(cif::category& cat);

	cif::datablock& mDb;
	std::streambuf* mIs;
	Token mLookahead;
	std::stack<unicode> mBuffer;
	std::string mToken;
	double mTokenFloat;
	int64_t mTokenInteger;
};

// -----------------------------------------------------------------------

uint8_t Parser::GetNextByte()
{
	int result = mIs->sbumpc();

	if (result == std::streambuf::traits_type::eof())
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
				throw std::runtime_error("Invalid utf-8");
			result = ((result & 0x01F) << 6) | (ch[0] & 0x03F);
		}
		else if ((result & 0x0F0) == 0x0E0)
		{
			ch[0] = GetNextByte();
			ch[1] = GetNextByte();
			if ((ch[0] & 0x0c0) != 0x080 or (ch[1] & 0x0c0) != 0x080)
				throw std::runtime_error("Invalid utf-8");
			result = ((result & 0x00F) << 12) | ((ch[0] & 0x03F) << 6) | (ch[1] & 0x03F);
		}
		else if ((result & 0x0F8) == 0x0F0)
		{
			ch[0] = GetNextByte();
			ch[1] = GetNextByte();
			ch[2] = GetNextByte();
			if ((ch[0] & 0x0c0) != 0x080 or (ch[1] & 0x0c0) != 0x080 or (ch[2] & 0x0c0) != 0x080)
				throw std::runtime_error("Invalid utf-8");
			result = ((result & 0x007) << 18) | ((ch[0] & 0x03F) << 12) | ((ch[1] & 0x03F) << 6) | (ch[2] & 0x03F);

			if (result > 0x10ffff)
				throw std::runtime_error("invalid utf-8 character (out of range)");
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
				throw std::runtime_error("character " + zeep::to_hex(result) + " is not allowed");

			// surrogate support
			else if (result >= 0x0D800 and result <= 0x0DBFF)
			{
				unicode uc2 = GetNextChar();
				if (uc2 >= 0x0DC00 and uc2 <= 0x0DFFF)
					result = (result - 0x0D800) * 0x400 + (uc2 - 0x0DC00) + 0x010000;
				else
					throw std::runtime_error("leading surrogate character without trailing surrogate character");
			}
			else if (result >= 0x0DC00 and result <= 0x0DFFF)
				throw std::runtime_error("trailing surrogate character without a leading surrogate");
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
	mBuffer.push(zeep::pop_last_char(mToken));
}

bool is_name_start_char(unicode uc)
{
	return
		(uc >= L'A' and uc <= L'Z') or
		uc == L'_' or
		(uc >= L'a' and uc <= L'z') or
		(uc >= 0x0C0 and uc <= 0x0D6) or
		(uc >= 0x0D8 and uc <= 0x0F6) or
		(uc >= 0x0F8 and uc <= 0x02FF) or
		(uc >= 0x0370 and uc <= 0x037D) or
		(uc >= 0x037F and uc <= 0x01FFF) or
		(uc >= 0x0200C and uc <= 0x0200D) or
		(uc >= 0x02070 and uc <= 0x0218F) or
		(uc >= 0x02C00 and uc <= 0x02FEF) or
		(uc >= 0x03001 and uc <= 0x0D7FF) or
		(uc >= 0x0F900 and uc <= 0x0FDCF) or
		(uc >= 0x0FDF0 and uc <= 0x0FFFD) or
		(uc >= 0x010000 and uc <= 0x0EFFFF);	
}

bool is_name_char(unicode uc)
{
	return
		(uc >= '0' and uc <= '9') or
		uc == 0x0B7 or
		is_name_start_char(uc) or
		(uc >= 0x00300 and uc <= 0x0036F) or
		(uc >= 0x0203F and uc <= 0x02040);
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
		EscapeHex4,

		Less,
		Greater
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
			case '=':
				token = Token::eq_;
				break;
			case '<':
				state = State::Less;
				break;
			case '>':
				state = State::Greater;
				break;
			case ' ':
			case '\n':
			case '\r':
			case '\t':
				mToken.clear();
				break;
			case '\'':
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
				else if (is_name_start_char(ch))
					state = State::Literal;
				else
					throw std::runtime_error("invalid character (" + zeep::to_hex(ch) + "/'" + (isprint(ch) ? static_cast<char>(ch) : '.') + "') in command");
			}
			break;
		
		case State::Less:
			if (ch == '=')
				token = Token::le_;
			else if (ch == '>')
				token = Token::ne_;
			else
			{
				Retract();
				token = Token::lt_;
			}
			break;

		case State::Greater:
			if (ch == '=')
				token = Token::ge_;
			else
			{
				Retract();
				token = Token::gt_;
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
				throw std::runtime_error("invalid character '-' in command");
			break;

		case State::NegativeZero:
			if (ch >= '0' or ch <= '9')
				throw std::runtime_error("invalid number in command, should not start with zero");
			token = Token::number;
			break;

		case State::Zero:
			if (ch >= '0' or ch <= '9')
				throw std::runtime_error("invalid number in command, should not start with zero");
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
				throw std::runtime_error("invalid floating point format in command");
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
			if (not is_name_char(ch))
			{
				Retract();

					 if (iequals(mToken, "SELECT")) 	token = Token::select;
				else if (iequals(mToken, "DISTINCT")) 	token = Token::distinct;
				else if (iequals(mToken, "FROM")) 		token = Token::from;
				else if (iequals(mToken, "UPDATE")) 	token = Token::update;
				else if (iequals(mToken, "SET")) 		token = Token::set;
				else if (iequals(mToken, "WHERE")) 		token = Token::where;
				else if (iequals(mToken, "AND")) 		token = Token::and_;
				else if (iequals(mToken, "OR")) 		token = Token::or_;
				else if (iequals(mToken, "NOT")) 		token = Token::not_;
				else if (iequals(mToken, "INSERT")) 	token = Token::insert;
				else if (iequals(mToken, "DELETE")) 	token = Token::delete_;
				else if (iequals(mToken, "INTO")) 		token = Token::into;
				else if (iequals(mToken, "VALUES")) 	token = Token::values;
				else if (iequals(mToken, "IS")) 		token = Token::is_;
				else if (iequals(mToken, "NULL")) 		token = Token::null_;
				else									token = Token::ident;
			}
			break;
		
		case State::String:
			if (ch == '\'')
			{
				token = Token::string;
				mToken.pop_back();
			}
			else if (ch == 0)
				throw std::runtime_error("Invalid unterminated std::string in command");
			else if (ch == '\\')
			{
				state = State::Escape;
				mToken.pop_back();
			}
			break;
		
		case State::Escape:
			switch (ch)
			{
				case '\'':	
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
					throw std::runtime_error("Invalid escape sequence in command (\\" + std::string{static_cast<char>(ch)} + ')');
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
				throw std::runtime_error("Invalid hex sequence in command");
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
				throw std::runtime_error("Invalid hex sequence in command");
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
				throw std::runtime_error("Invalid hex sequence in command");
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
				throw std::runtime_error("Invalid hex sequence in command");
			mToken.pop_back();
			zeep::append(mToken, hx);
			state = State::String;
			break;
		}
	}

	return token;
}

void Parser::Match(Token expected)
{
	if (mLookahead != expected)
		throw std::runtime_error("Syntax error in command, expected " + Describe(expected) + " but found " + Describe(mLookahead) + " (" + mToken + ")");
	
	mLookahead = GetNextToken();
}

StatementPtr Parser::Parse(std::streambuf* is)
{
	mIs = is;

	mLookahead = GetNextToken();
	std::shared_ptr<StatementList> result(new StatementList());

	while (mLookahead != Token::eoln)
	{
		auto stmt = ParseStatement();
		result->Add(stmt);
	}
	
	return result;
}

// -----------------------------------------------------------------------

StatementPtr Parser::ParseStatement()
{
	StatementPtr result;

	switch (mLookahead)
	{
		case Token::select:
			Match(Token::select);
			result = ParseSelect();
			break;
		
		case Token::delete_:
			Match(Token::delete_);
			result = ParseDelete();
			break;
		
		case Token::update:
			Match(Token::update);
			result = ParseUpdate();
			break;
		
		default:
			// force error
			Match(Token::select);
	}
	
	Match(Token::semicolon);

	return result;
}

// -----------------------------------------------------------------------

StatementPtr Parser::ParseSelect()
{
	bool distinct = false;
	if (mLookahead == Token::distinct)
	{
		distinct = true;
		Match(Token::distinct);
	}

	auto items = ParseItemList();

	Match(Token::from);

	std::string cat = mToken;
	Match(Token::ident);

	auto category = mDb.get(cat);
	if (category == nullptr)
		throw std::runtime_error("Category " + cat + " is not defined in this file");

	auto cv = category->get_cat_validator();
	if (cv != nullptr)
	{
		std::vector<std::string> nItems;

		for (auto item: items)
		{
			if (item == "*")
				transform(cv->m_item_validators.begin(), cv->m_item_validators.end(), back_inserter(nItems),
					[cat](auto iv) { return iv.m_tag; });
			else
				nItems.push_back(item);
		}

		swap(items, nItems);

		items.erase(remove_if(items.begin(), items.end(), [category](auto item) { return not category->has_column(item); }), items.end());

		for (auto item: items)
		{
			auto iv = cv->get_validator_for_item(item);
			if (iv == nullptr)
				throw std::runtime_error("Item " + item + " is not defined in the PDBx dictionary for category " + cat);
		}	
	}

	if (mLookahead == Token::where)
	{
		Match(Token::where);
		return StatementPtr{ new SelectStatement(*category, distinct, std::move(items), ParseNotWhereClause(*category)) };
	}
	else
		return StatementPtr{ new SelectStatement(*category, distinct, std::move(items), cif::all()) };
}

// -----------------------------------------------------------------------

StatementPtr Parser::ParseDelete()
{
	Match(Token::from);

	std::string cat = mToken;
	Match(Token::ident);

	auto category = mDb.get(cat);
	if (category == nullptr)
		throw std::runtime_error("Category " + cat + " is not defined in this file");

	if (mLookahead == Token::where)
	{
		Match(Token::where);
		return StatementPtr{ new DeleteStatement(*category, ParseNotWhereClause(*category)) };
	}
	else
		return StatementPtr{ new DeleteStatement(*category, cif::all()) };
}

// -----------------------------------------------------------------------

StatementPtr Parser::ParseUpdate()
{
	std::string cat = mToken;
	Match(Token::ident);

	auto category = mDb.get(cat);
	if (category == nullptr)
		throw std::runtime_error("Category " + cat + " is not defined in this file");

	auto cv = category->get_cat_validator();

	Match(Token::set);

	std::vector<std::pair<std::string,std::string>> itemValuePairs;
	for (;;)
	{
		std::string item = mToken;
		Match(Token::ident);

		auto iv = cv ? cv->get_validator_for_item(item) : nullptr;
		if (cv and iv == nullptr)
			throw std::runtime_error("Invalid item '" + item + "' for category '" + cat + '\'');
		
		Match(Token::eq_);

		std::string value = mToken;
		switch (mLookahead)
		{
			case Token::integer:
			case Token::number:
			case Token::string:
				Match(mLookahead);
				break;
			default:
				Match(Token::string);
		}

		if (iv)
			iv->operator()(value);

		itemValuePairs.emplace_back(item, value);

		if (mLookahead == Token::comma)
		{
			Match(Token::comma);
			continue;
		}

		break;
	}

	if (mLookahead == Token::where)
	{
		Match(Token::where);
		return StatementPtr{ new UpdateStatement(*category, std::move(itemValuePairs), ParseNotWhereClause(*category)) };
	}
	else
		return StatementPtr{ new UpdateStatement(*category, std::move(itemValuePairs), cif::all()) };
}

// -----------------------------------------------------------------------

std::vector<std::string> Parser::ParseItemList()
{
	std::vector<std::string> items;

	for (;;)
	{
		if (mLookahead == Token::asterisk)
		{
			Match(Token::asterisk);
			items.push_back("*");
		}
		else
		{
			items.push_back(mToken);
			Match(Token::ident);
		}

		if (mLookahead == Token::comma)
		{
			Match(Token::comma);
			continue;
		}
		
		break;
	}

	return items;
}

// -----------------------------------------------------------------------

cif::condition Parser::ParseNotWhereClause(cif::category& cat)
{
	cif::condition result;

	if (mLookahead == Token::not_)
	{
		Match(Token::not_);
		result = not ParseNotWhereClause(cat);
	}
	else if (mLookahead == Token::braceopen)
	{
		Match(Token::braceopen);
		result = ParseNotWhereClause(cat);
		Match(Token::braceclose);
	}
	else
	{
		result = ParseWhereClause(cat);

		for (;;)
		{
			if (mLookahead == Token::and_)
			{
				Match(Token::and_);
				result = std::move(result) and ParseNotWhereClause(cat);
				continue;
			}

			if (mLookahead == Token::or_)
			{
				Match(Token::or_);
				result = std::move(result) or ParseNotWhereClause(cat);
				continue;
			}

			break;
		}
	}
	return result;
}

// -----------------------------------------------------------------------

cif::condition Parser::ParseWhereClause(cif::category& cat)
{
	std::string item = mToken;
	Match(Token::ident);

	auto cv = cat.get_cat_validator();
	if (cv != nullptr and cv->get_validator_for_item(item) == nullptr)
	{
		throw std::runtime_error("Invalid item '" + item + "' for category '" + cat.name() + "' in where clause");
	}

	if (mLookahead == Token::is_)
	{
		Match(mLookahead);

		if (mLookahead == Token::not_)
		{
			Match(mLookahead);
			Match(Token::null_);
			return cif::key(item) != cif::null;
		}
		else
		{
			Match(Token::null_);
			return cif::key(item) == cif::null;
		}
	}
	else
	{
		if (mLookahead < Token::eq_ or mLookahead > Token::ne_)
			Match(Token::eq_);
		
		auto oper = mLookahead;
		Match(mLookahead);
		
		cif::condition c;
		std::string value = mToken;

		switch (mLookahead)
		{
			case Token::integer:
			case Token::number:
			case Token::string:
				Match(mLookahead);
				break;
			default:
				Match(Token::string);
		}

		switch (oper)
		{
			case Token::eq_:	return cif::key(item) == value;
			case Token::lt_:	return cif::key(item) <  value;
			case Token::le_:	return cif::key(item) <= value;
			case Token::gt_:	return cif::key(item) >  value;
			case Token::ge_:	return cif::key(item) >= value;
			case Token::ne_:	return cif::key(item) != value;
			default:			throw std::logic_error("should never happen");
		}
	}
}

}

// -----------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("mmCQL [options] input output");
	visible_options.add_options()
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")

		("force",										"Force writing of output file, even if it is the same as the input file")

		("script,f",     po::value<std::string>(),   		"Read commands from script");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",		po::value<std::string>(),			"Input file")
		("output,o",	po::value<std::string>(),			"Output file")
		("debug,d",		po::value<int>(),				"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		write_version_string(std::cout, vm.count("verbose"));
		exit(0);
	}

	if (vm.count("help") or not vm.count("input"))
	{
		std::cerr << visible_options << std::endl;
		exit(vm.count("help") != 0);
	}

	if (vm.count("output") and vm["output"].as<std::string>() == vm["input"].as<std::string>() and vm.count("force") == 0)
	{
		std::cerr << "Cowardly refusing to overwrite input file (specify --force to force overwriting)" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	auto input = vm["input"].as<std::string>();
	cif::file file{fs::path(input)};

	cql::Parser parser(file.front());

	if (vm.count("script"))
	{
		std::ifstream cmdFile(vm["script"].as<std::string>());
		if (not cmdFile.is_open())
			throw std::runtime_error("Failed to open command file " + vm["script"].as<std::string>());

		auto stmt = parser.Parse(cmdFile.rdbuf());
		if (stmt)
			stmt->Execute();
	}
	else
	{
		std::string cmd;
		while (std::getline(std::cin, cmd))
		{
			try
			{
				std::istringstream is(cmd);

				auto stmt = parser.Parse(is.rdbuf());
				if (stmt)
					stmt->Execute();
			}
			catch(const std::exception& e)
			{
				std::cerr << e.what() << std::endl;
			}
		}
	}

	if (vm.count("output"))
		file.save(vm["output"].as<std::string>());

	return 0;	
}

