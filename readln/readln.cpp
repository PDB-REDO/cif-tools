/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 Maarten L. Hekkelman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "readln.hpp"

#include <cassert>
#include <cerrno>
#include <deque>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <csignal>
#include <cstdlib>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <termios.h>
#include <unistd.h>

// --------------------------------------------------------------------

namespace readln
{

class editor
{
	enum ControlCharacter : uint8_t // NOLINT
	{
		// clang-format off
		NUL, SOH, STX, ETX, EOT, ENQ, ACK, BEL,  BS,  HT,  LF,  VT,  FF,  CR,  SO,  SI,
		DLE, DC1, DC2, DC3, DC4, NAK, SYN, ETB, CAN,  EM, SUB, ESC,  FS,  GS,  RS,  US,
		SP = 0x20,
		DEL = 0x7f,
		IND = 0x84, NEL, SSA, ESA, HTS, HTJ, VTS, PLD, PLU, RI, SS2, SS3,
		DCS, PU1, PU2, STS, CCH, MW, SPA, EPA, CSI, ST, OSC, PM, APC
		// clang-format on
	};

  public:
	editor(std::string_view prompt, int in_fd, int out_fd);

	~editor();

	bool getline(std::string &line);

	static void load_history(const std::filesystem::path &file);
	static void save_history(const std::filesystem::path &file);
	static void set_history_max_len(size_t len);
	static void add_history(std::string_view s);

  private:
	enum class action
	{
		CONTINUE,
		DONE,
		ABORT
	};

	uint8_t read_char();

	enum command : int
	{
		IGNORE = 256,

		ABORT,
		ENTER,
		INSERT,
		BACKSPACE,
		DELETE,

		CURSOR_UP,
		CURSOR_DOWN,
		CURSOR_RIGHT,
		CURSOR_LEFT,
		CURSOR_END,
		CURSOR_HOME,

		CURSOR_WORD_RIGHT,
		CURSOR_WORD_LEFT,

		// SELECT_CHAR_RIGHT,
		// SELECT_CHAR_LEFT,
		// SELECT_WORD_RIGHT,
		// SELECT_WORD_LEFT,

		ERASE_TO_EOLN,
		ERASE_TO_BOLN,
		ERASE_LINE,
		ERASE_WORD_LEFT,

		CLEAR_SCREEN

	};

	std::tuple<command, char> read_command();

	void repaint();
	void place_cursor();
	void beep();

	int word_left(int pos);
	int word_right(int pos);

	static constexpr const int kWordStateTable[2][2] = {
		{ 0, 1 },
		{ -1, 1 }
	};

	std::string m_prompt;
	int m_in_fd, m_out_fd;
	termios m_termios_state;
	bool m_raw_mode = false;

	std::string m_buffer;
	int m_width, m_cursor, m_offset;

	void win_resized();
	static void win_resized_handler(int signo, siginfo_t *info, void *context);
	struct sigaction m_old_sigaction;
	bool m_sig_action_set = false;
	static bool s_signaled;

	// history
	static std::deque<std::string> s_history;
	static std::deque<std::string>::iterator s_history_ptr;
	static size_t s_max_history_len;
};

// --------------------------------------------------------------------

bool editor::s_signaled = false;

std::deque<std::string> editor::s_history; // NOLINT
std::deque<std::string>::iterator editor::s_history_ptr;
size_t editor::s_max_history_len = 100;

// --------------------------------------------------------------------

editor::editor(std::string_view prompt, int in_fd, int out_fd)
	: m_prompt(prompt)
	, m_in_fd(in_fd)
	, m_out_fd(out_fd)
{
	auto r = tcgetattr(m_in_fd, &m_termios_state);
	if (r != 0)
		throw std::system_error(errno, std::system_category(), "Not a terminal?");

	auto new_state = m_termios_state;
	new_state.c_iflag &=
		~(IGNBRK | BRKINT | PARMRK | ISTRIP | INLCR | IGNCR | ICRNL | IXON);
	new_state.c_oflag &= ~OPOST;
	new_state.c_lflag &= ~(ECHO | ECHONL | ICANON | ISIG | IEXTEN);
	new_state.c_cflag &= ~(CSIZE | PARENB);
	new_state.c_cflag |= CS8;

	new_state.c_cc[VMIN] = 1;
	new_state.c_cc[VTIME] = 1;

	r = tcsetattr(m_in_fd, TCSAFLUSH, &new_state);
	if (r != 0)
		throw std::system_error(errno, std::system_category(), "Not a terminal?");

	m_raw_mode = true;

	// --------------------------------------------------------------------

	struct winsize w;
	::ioctl(m_in_fd, TIOCGWINSZ, &w);
	m_width = w.ws_col - m_prompt.length();
	m_cursor = 0;
	m_offset = 0;

	struct sigaction act = {};

	act.sa_flags = SA_SIGINFO /*  | SA_UNSUPPORTED | SA_EXPOSE_TAGBITS */;
	act.sa_sigaction = &win_resized_handler;
	if (sigaction(SIGWINCH, &act, &m_old_sigaction) == 0)
		m_sig_action_set = true;
}

editor::~editor()
{
	if (m_sig_action_set)
		sigaction(SIGWINCH, &m_old_sigaction, nullptr);

	if (m_raw_mode)
		tcsetattr(m_in_fd, TCSAFLUSH, &m_termios_state);
}

void editor::load_history(const std::filesystem::path &file)
{
	s_history.clear();

	std::ifstream f(file);
	if (f.is_open())
	{
		std::string line;
		while (std::getline(f, line))
		{
			s_history.emplace_back(line);
			if (s_history.size() > s_max_history_len)
				s_history.pop_front();
		}
	}
}

void editor::save_history(const std::filesystem::path &file)
{
	std::ofstream f(file);
	if (f.is_open())
	{
		for (auto &line : s_history)
			f << line << '\n';
	}
}

void editor::set_history_max_len(size_t len)
{
	s_max_history_len = len;
	if (s_history.size() > s_max_history_len)
		s_history.erase(s_history.begin(), s_history.begin() + s_history.size() - s_max_history_len);
}

void editor::add_history(std::string_view s)
{
	if (s_history.empty() or s_history.back() != s)
		s_history.emplace_back(s);
}

void editor::win_resized_handler(int signo, siginfo_t *info, void *context)
{
	s_signaled = true;
}

uint8_t editor::read_char()
{
	char result;
	for (;;)
	{
		auto r = read(m_in_fd, &result, 1);
		if (r < 0)
		{
			if (errno == EINTR)
			{
				if (s_signaled)
				{
					struct winsize w;
					::ioctl(m_in_fd, TIOCGWINSZ, &w);
					m_width = w.ws_col - m_prompt.length();
					s_signaled = false;

					repaint();
				}

				continue;
			}

			throw std::system_error(errno, std::system_category(), "reading character from terminal");
		}
		break;
	}

	return result;
}

auto editor::read_command() -> std::tuple<command, char>
{
	command result = IGNORE;
	char ch = 0;

	enum class State
	{
		START,
		ESCAPE,
		CSI,
		CAN,
		DCS_APC_OSC,
		CharSet_TransM
	} state = State::START;

	std::vector<int> args;
	int csiCmd = 0;

	while (result == IGNORE)
	{
		auto c = read_char();

		switch (state)
		{
			case State::START:
				switch (c)
				{
					case ETX:
						result = ABORT;
						break;

					case EOT:
						if (m_buffer.empty())
							result = ABORT;
						else
							beep();
						break;

					case CR:
					case LF:
						result = ENTER;
						break;

					case ESC:
						state = State::ESCAPE;
						break;

					case BEL:
						beep();
						break;

					case FF:
						result = CLEAR_SCREEN;
						break;

					case DEL:
						result = BACKSPACE;
						break;

					case VT:
						result = ERASE_TO_EOLN;
						break;

					case NAK:
						result = ERASE_LINE;
						break;

					case CAN:
						state = State::CAN;
						break;

					case ETB:
						result = ERASE_WORD_LEFT;
						break;

					default:
						if ((c >= ' ' and c <= '/') or (c >= '0' and c <= '9') or (c >= ':' and c <= '~'))
						{
							result = INSERT;
							ch = static_cast<char>(c);
						}
						else
							result = static_cast<command>(c);
						break;
				}

				break;

			case State::CAN:
				if (c == DEL)
					result = ERASE_TO_BOLN;
				else
					state = State::START;
				break;

			case State::ESCAPE:
				switch (c)
				{
					case '[':
						state = State::CSI;
						csiCmd = 0;
						args = { 0 };
						break;

					case '_':
						state = State::DCS_APC_OSC;
						break;

					case '*':
					case '+':
					case ' ':
						state = State::CharSet_TransM;
						break;

					default:
						state = State::START;
				}
				break;

			case State::CSI:
				if (c >= '0' and c <= '9')
					args.back() = args.back() * 10 + c - '0';
				else if (c == ';')
					args.push_back(0);
				else if (c == ' ' or c == '?')
					csiCmd = csiCmd << 8 | c;
				else
				{
					state = State::START;
					if (c < '0' or c > '~')
						state = State::START;
					else
					{
						csiCmd = csiCmd << 8 | c;
						switch (csiCmd)
						{
							case 'A':
								result = CURSOR_UP;
								break;

							case 'B':
								result = CURSOR_DOWN;
								break;

							case 'C':
								switch (args[1])
								{
									// case 2:
									// 	result = SELECT_CHAR_RIGHT;
									// 	break;
									case 5:
										result = CURSOR_WORD_RIGHT;
										break;
									// case 6:
									// 	result = SELECT_WORD_RIGHT;
									// 	break;
									default:
										result = CURSOR_RIGHT;
										break;
								}
								break;

							case 'D':
								switch (args[1])
								{
									// case 2:
									// 	result = SELECT_CHAR_LEFT;
									// 	break;
									case 5:
										result = CURSOR_WORD_LEFT;
										break;
									// case 6:
									// 	result = SELECT_WORD_LEFT;
									// 	break;
									default:
										result = CURSOR_LEFT;
										break;
								}
								break;

							case 'F':
								result = CURSOR_END;
								break;

							case 'H':
								result = CURSOR_HOME;
								break;

							case 'J':
								result = ERASE_TO_EOLN;
								break;

							case 'K':
								result = ERASE_TO_BOLN;
								break;

							case 'Y':
								read_char();
								read_char();
								state = State::START;
								break;

							case '~':
								switch (args[0])
								{
									case 1:
										result = CURSOR_HOME;
										break;

									case 3:
										result = DELETE;
										break;

									case 4:
										result = CURSOR_END;
										break;
									
									default:
										break;
								}
								break;

							default:
								break;
						}
					}
				}
				break;

			case State::DCS_APC_OSC:
				if (c == ST or c == BEL)
					state = State::START;
				break;

			case State::CharSet_TransM:
				state = State::START;
				break;
		}
	}

	return { result, ch };
}

bool editor::getline(std::string &line)
{
	line.clear();
	m_buffer.clear();
	m_cursor = 0;
	m_offset = 0;
	s_history.emplace_back("");
	s_history_ptr = s_history.end() - 1;

	repaint();

	bool result = false;

	for (;;)
	{
		bool dirty = false;
		auto oldCursor = m_cursor;
		auto oldOffset = m_offset;

		action act = action::CONTINUE;

		switch (const auto &[cmd, ch] = read_command(); cmd)
		{
			case IGNORE:
				break;

			case ABORT:
				act = action::ABORT;
				break;

			case ENTER:
				act = action::DONE;
				break;

			case INSERT:
				m_buffer.insert(m_buffer.begin() + m_cursor, ch);
				m_cursor += 1;
				dirty = true;
				break;

			case BACKSPACE:
				if (m_cursor > 0 and std::cmp_less_equal(m_cursor, m_buffer.length()))
				{
					m_buffer.erase(m_buffer.begin() + m_cursor - 1);
					--m_cursor;
					dirty = true;
				}
				break;

			case DELETE:
				if (m_cursor >= 0 and std::cmp_less(m_cursor + 1, m_buffer.length()))
				{
					m_buffer.erase(m_buffer.begin() + m_cursor);
					dirty = true;
				}
				break;

			case CURSOR_LEFT:
				if (m_cursor > 0)
					--m_cursor;
				break;

			case CURSOR_RIGHT:
				if (std::cmp_less(m_cursor, m_buffer.length()))
					++m_cursor;
				break;

			case CURSOR_HOME:
				m_cursor = 0;
				break;

			case CURSOR_END:
				m_cursor = m_buffer.length();
				break;

			case CURSOR_WORD_RIGHT:
				m_cursor = word_right(m_cursor);
				break;

			case CURSOR_WORD_LEFT:
				m_cursor = word_left(m_cursor);
				break;

				// case SELECT_CHAR_RIGHT:
				// 	break;

				// case SELECT_CHAR_LEFT:
				// 	break;

				// case SELECT_WORD_RIGHT:
				// 	break;

				// case SELECT_WORD_LEFT:
				// 	break;

			case ERASE_WORD_LEFT:
			{
				auto p = word_left(m_cursor);
				m_buffer.erase(p, m_cursor - p);
				dirty = true;
				m_cursor = p;
				break;
			}

			case ERASE_TO_EOLN:
				if (std::cmp_less(m_cursor, m_buffer.length()))
				{
					m_buffer.erase(m_cursor, std::string::npos);
					dirty = true;
				}
				break;

			case ERASE_TO_BOLN:
				if (m_cursor > 0)
				{
					m_buffer.erase(0, m_cursor);
					m_cursor = 0;
					dirty = true;
				}
				break;

			case ERASE_LINE:
				m_buffer.clear();
				m_cursor = 0;
				dirty = true;
				break;

			case CLEAR_SCREEN:
				write(m_out_fd, "\x1b[2J\x1b[H", 7);
				m_cursor = 0;
				dirty = true;
				break;

			case CURSOR_UP:
				if (s_history_ptr != s_history.begin())
				{
					std::swap(m_buffer, *s_history_ptr);
					--s_history_ptr;
					m_buffer = *s_history_ptr;
					m_cursor = m_buffer.length();
					dirty = true;
				}
				break;

			case CURSOR_DOWN:
				if (s_history_ptr + 1 < s_history.end())
				{
					std::swap(m_buffer, *s_history_ptr);
					++s_history_ptr;
					m_buffer = *s_history_ptr;
					m_cursor = m_buffer.length();
					dirty = true;
				}
				break;

			default:
				beep();
				break;
		}

		if (m_offset < m_cursor - m_width + 1)
			m_offset = m_cursor - m_width + 1;
		else if (m_offset > m_cursor)
			m_offset = m_cursor;

		if (dirty or m_offset != oldOffset)
			repaint();
		else if (m_cursor != oldCursor)
			place_cursor();

		if (act == action::CONTINUE)
			continue;

		if (act == action::DONE)
		{
			line = m_buffer;
			result = true;
		}

		break;
	}

	s_history.pop_back();

	write(m_out_fd, "\r\n", 2);

	return result;
}

void editor::repaint()
{
	// erase line
	write(m_out_fd, "\x1b[2K", 4);
	// move to start
	write(m_out_fd, "\x1b[1G", 4);
	// write prompt
	write(m_out_fd, m_prompt.data(), m_prompt.length());
	// write out buffer, but only the part that fits

	auto o = m_offset;
	assert(o >= 0 and o <= m_buffer.length());
	auto l = m_buffer.length() - m_offset;

	if (std::cmp_greater(l, m_width))
		l = m_width;

	write(m_out_fd, m_buffer.data() + o, l);

	place_cursor();
}

void editor::place_cursor()
{
	// set cursor

	auto pos = m_prompt.length() + m_cursor - m_offset;
	auto pcmd = std::format("\x1b[{}G", pos + 1);
	write(m_out_fd, pcmd.data(), pcmd.length());
}

void editor::beep()
{
	write(m_out_fd, "\x07", 1);
}

int editor::word_left(int pos)
{
	auto c = m_buffer.begin() + pos;
	int8_t state = 0;

	while (state >= 0)
	{
		if (c == m_buffer.begin())
			return 0;

		uint8_t ch = *--c;
		int wc = isalnum(ch) or ch == '_' ? 1 : 0;
		state = kWordStateTable[state][wc];
	}

	return c - m_buffer.begin() + 1;
}

int editor::word_right(int pos)
{
	auto c = m_buffer.begin() + pos;
	int8_t state = 0;

	while (state >= 0)
	{
		if (c == m_buffer.end())
			return m_buffer.length();

		uint8_t ch = *c++;
		int wc = isalnum(ch) or ch == '_' ? 1 : 0;
		state = kWordStateTable[state][wc];
	}

	return c - m_buffer.begin() - 1;
}

// --------------------------------------------------------------------

bool getline(std::string prompt, std::string &line)
{
	std::cout.flush();

	editor editor(prompt, STDIN_FILENO, STDOUT_FILENO);
	bool result = editor.getline(line);
	return result;
}

void load_history(const std::filesystem::path &file)
{
	editor::load_history(file);
}

void save_history(const std::filesystem::path &file)
{
	editor::save_history(file);
}

void set_history_max_len(size_t len)
{
	editor::set_history_max_len(len);
}

void add_history(std::string_view s)
{
	editor::add_history(s);
}

} // namespace readln