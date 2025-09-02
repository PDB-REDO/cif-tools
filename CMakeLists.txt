# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAcGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.23)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# set the project name
project(cif-tools VERSION 1.0.12 LANGUAGES CXX)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CPM)
include(GNUInstallDirs)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

if(NOT "$ENV{CCP4}" STREQUAL "")
	set(CCP4 $ENV{CCP4})
	list(PREPEND CMAKE_MODULE_PATH "${CCP4}/Lib")
	list(APPEND CMAKE_PREFIX_PATH ${CCP4})

	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		set(CMAKE_PREFIX_PATH ${CCP4})
	endif()
endif()

if(WIN32)
	if(${CMAKE_SYSTEM_VERSION} GREATER_EQUAL 10) # Windows 10
		add_definitions(-D _WIN32_WINNT=0x0A00)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.3) # Windows 8.1
		add_definitions(-D _WIN32_WINNT=0x0603)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.2) # Windows 8
		add_definitions(-D _WIN32_WINNT=0x0602)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.1) # Windows 7
		add_definitions(-D _WIN32_WINNT=0x0601)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.0) # Windows Vista
		add_definitions(-D _WIN32_WINNT=0x0600)
	else() # Windows XP (5.1)
		add_definitions(-D _WIN32_WINNT=0x0501)
	endif()

	# Man, this is 2024 we're living in...
	add_definitions(-DNOMINMAX)

	# We do not want to write an export file for all our symbols...
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive- /bigobj)
	add_link_options(/NODEFAULTLIB:library)

	# This is dubious...
	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("${PROJECT_SOURCE_DIR}/src")

# Optionally use mrc to create resources
find_package(mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "cif-tools: Not using resources since mrc was not found")
	set(USE_RSRC FALSE)
endif()

if(USE_RSRC)
	add_compile_definitions(USE_RSRC)
endif()

if(NOT(mcfp_FOUND OR libmcfp_FOUND OR TARGET mcfp))
	set(CPM_USE_LOCAL_PACKAGES ON)
	CPMAddPackage("https://forge.hekkelman.net/maarten/mcfp.git@1.4.2")
endif()

# find_package(cifpp 9 QUIET)

if(NOT(cifpp_FOUND OR TARGET cifpp))
	add_subdirectory(libcifpp EXCLUDE_FROM_ALL)
	set(cifpp_FOUND TRUE)
	set(CIFPP_SHARE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libcifpp/rsrc)
endif()

list(APPEND programs
	pdb2cif
	cif2pdb
	cif-diff
	cif-grep
	cif-merge
	cif-validate
	mmCQL
)

add_library(prmain OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/src/pr-main.cpp")
target_link_libraries(prmain cifpp::cifpp)

foreach(PROGRAM IN LISTS programs)
	add_executable(${PROGRAM} ${PROJECT_SOURCE_DIR}/src/${PROGRAM}.cpp)

	if(USE_RSRC)
		mrc_target_resources(${PROGRAM}
			${CIFPP_SHARE_DIR}/mmcif_ma.dic
			${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
			${CIFPP_SHARE_DIR}/mmcif_ddl.dic)
	endif()

	target_link_libraries(${PROGRAM} PRIVATE cifpp::cifpp mcfp::mcfp $<TARGET_OBJECTS:prmain>)

	install(TARGETS ${PROGRAM} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

	install(FILES doc/${PROGRAM}.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endforeach()
