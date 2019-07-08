/* 
   Created by: Maarten L. Hekkelman
   Date: woensdag 21 september, 2016
*/

#include "libpr.h"

#include <iostream>

#include <boost/filesystem.hpp>

#include "centrifuge.h"
#include "cif-diff.h"
#include "cif-grep.h"
#include "cif-test.h"
#include "cif-merge.h"
#include "stripper.h"
#include "extractor.h"
#include "tls-interpreter.h"

using namespace std;
namespace fs = boost::filesystem;

const string kAppName = "pr-driver";

// --------------------------------------------------------------------
// global variables

int VERBOSE;

// --------------------------------------------------------------------
//    main

int main(int argc, char* argv[])
{
	int result = 0;
	
	try
	{
		string command = fs::path(argv[0]).leaf().string();

		if (command == kAppName and argc > 1)
		{
			command = argv[1];
			++argv;
			--argc;
		}
		
		if (command == "centrifuge-learn")
			result = centrifuge_learn(argc, argv);
		else if (command == "centrifuge-predict")
			result = centrifuge_predict(argc, argv);
		else if (command == "centrifuge-test")
			result = centrifuge_test(argc, argv);
		else if (command == "stripper")
			result = stripper(argc, argv);
		else if (command == "extractor")
			result = extractor(argc, argv);
		else if (command == "edia-test")
			result = edia_test(argc, argv);
		else if (command == "tls-interpreter")
			result = tls_interpreter_main(argc, argv);
		else if (command == "cif-diff")
			result = cif_diff(argc, argv);
		else if (command == "cif-grep")
			result = cif_grep(argc, argv);
		else if (command == "cif-drop")
			result = cif_drop(argc, argv);
		else if (command == "cif-merge")
			result = cif_merge(argc, argv);
		else if (command == "cif-test")
			result = cif_test(argc, argv);
		else if (command == "pdb2cif")
			result = pdb2cif_test(argc, argv);
		else if (command == "cif2pdb")
			result = cif2pdb_test(argc, argv);
		else if (command == "pdb2cif-diff")
			result = cif_diff_test(argc, argv);
		else if (command == "cif2pdb-diff")
			result = pdb_diff_test(argc, argv);
		else if (command == "pdb2pdb-diff")
			result = pdb2pdb_diff_test(argc, argv);

// other commands here
		else
		{
			cout << "Usage: " << kAppName << " command [options]" << endl
				 << endl
				 << "Currently supported commands are:" << endl
				 << endl;
			
			for (auto cmd: { "centrifuge-learn", "centrifuge-predict", 
				"stripper",
				"cif-test", "cif-diff", "cif-merge", "cif-drop",
				"cif2pdb", "pdb2cif" })
			{
				cout << "\t" << cmd << endl;
			}

			exit(1);
		}
	}
	catch (const exception& e)
	{
		cerr << endl
			 << "pr exited with an exception: " << endl
			 << e.what() << endl;
		result = 1;	
	}
	
	catch (...)
	{
		cerr << endl
			 << "pr exited with an uncaught exception" << endl;
		result = 1;	
	}
	
	return result;
}
