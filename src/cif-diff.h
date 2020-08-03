#pragma once

#include <boost/filesystem/operations.hpp>
#include "mmcif/cif++.h"

int cif_diff(int argc, char* argv[]);
int cif_diff_test(int argc, char* argv[]);
int cif_drop(int argc, char* argv[]);

void compare_cifs_text(cif::file& a, cif::file& b, std::filesystem::path fa, bool icase, bool iwhite);
