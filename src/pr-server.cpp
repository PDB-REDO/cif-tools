/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 30 april, 2019
*/

#include "pdb-redo.h"

#include <csignal>
#include <pwd.h>
#include <sys/wait.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <thread>

#include <zeep/http/webapp.hpp>
#include <zeep/http/preforked-server.hpp>

#include "pr-server.hpp"

using namespace std;

namespace zh = zeep::http;

// --------------------------------------------------------------------
//
//	Daemonize
// 

void Daemonize()
{
	int pid = fork();
	
	if (pid == -1)
	{
		cerr << "Fork failed" << endl;
		exit(1);
	}
	
	if (pid != 0)
		_exit(0);

	if (setsid() < 0)
	{
		cerr << "Failed to create process group: " << strerror(errno) << endl;
		exit(1);
	}

	// it is dubious if this is needed:
	signal(SIGHUP, SIG_IGN);

	// fork again, to avoid being able to attach to a terminal device
	pid = fork();

	if (pid == -1)
		cerr << "Fork failed" << endl;

	if (pid != 0)
		_exit(0);

	// write our pid to the pid file
	ofstream pidFile(PID_FILE);
	pidFile << getpid() << endl;
	pidFile.close();

	if (chdir("/") != 0)
	{
		cerr << "Cannot chdir to /: " << strerror(errno) << endl;
		exit(1);
	}

	// close stdin
	close(STDIN_FILENO);
	open("/dev/null", O_RDONLY);
}

// --------------------------------------------------------------------
// 
//	OpenLogFile
//	

void OpenLogFile(const string& logFile)
{
	// open the log file
	int fd = open(logFile.c_str(), O_CREAT|O_APPEND|O_RDWR, 0644);
	if (fd < 0)
	{
		cerr << "Opening log file " << logFile << " failed" << endl;
		exit(1);
	}

	// redirect stdout and stderr to the log file
	dup2(fd, STDOUT_FILENO);
	dup2(fd, STDERR_FILENO);
}

// --------------------------------------------------------------------
// 
//	Main Loop
//

bool RunMainLoop(const string& address, uint16_t port, const string& user, CreateServer&& serverFactory)
{
	cout << "restarting server" << endl;

    sigset_t new_mask, old_mask;
    sigfillset(&new_mask);
    pthread_sigmask(SIG_BLOCK, &new_mask, &old_mask);

	zeep::http::preforked_server server([=]() -> zeep::http::server*
	{
		try
		{
			if (not user.empty())
			{
				struct passwd* pw = getpwnam(user.c_str());
				if (pw == NULL or setuid(pw->pw_uid) < 0)
				{
					cerr << "Failed to set uid to " << user << ": " << strerror(errno) << endl;
					exit(1);
				}
			}
			
			return serverFactory();
		}
		catch (const exception& e)
		{
			cerr << "Failed to launch server: " << e.what() << endl;
			exit(1);
		}
	});
	
	std::thread t(std::bind(&zeep::http::preforked_server::run, &server, address, port, 2));
	server.start();

	try
	{
		server.start();
	}
	catch (const exception& ex)
	{
		cerr << endl
			 << "Exception running server: " << endl
			 << ex.what() << endl
			 << endl;
	}

    pthread_sigmask(SIG_SETMASK, &old_mask, 0);

	// Wait for signal indicating time to shut down.
	sigset_t wait_mask;
	sigemptyset(&wait_mask);
	sigaddset(&wait_mask, SIGINT);
	sigaddset(&wait_mask, SIGHUP);
	sigaddset(&wait_mask, SIGQUIT);
	sigaddset(&wait_mask, SIGTERM);
	sigaddset(&wait_mask, SIGCHLD);
	pthread_sigmask(SIG_BLOCK, &wait_mask, 0);
	int sig = 0;
	sigwait(&wait_mask, &sig);

    pthread_sigmask(SIG_SETMASK, &old_mask, 0);

	server.stop();
	t.join();

	if (sig == SIGCHLD)
	{
		int status, pid;
		pid = waitpid(-1, &status, WUNTRACED);

		if (pid != -1 and WIFSIGNALED(status))
			cout << "child " << pid << " terminated by signal " << WTERMSIG(status) << endl;
	}

	return sig == SIGHUP;
}
