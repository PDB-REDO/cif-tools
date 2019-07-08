/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 30 april, 2019
*/

#pragma once

#include <zeep/http/server.hpp>

typedef zeep::http::server* (*CreateServer)();

// define this in your main code:
extern const char* PID_FILE;

void Daemonize();
void OpenLogFile(const std::string& logFile);
bool RunMainLoop(const std::string& address, uint16_t port,
	const std::string& user, CreateServer&& serverFactory);
