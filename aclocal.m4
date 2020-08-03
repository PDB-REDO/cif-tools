# generated automatically by aclocal 1.16.1 -*- Autoconf -*-

# Copyright (C) 1996-2018 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_boost_base.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_BASE([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for the Boost C++ libraries of a particular version (or newer)
#
#   If no path to the installed boost library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $BOOST_ROOT environment variable. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_CPPFLAGS) / AC_SUBST(BOOST_LDFLAGS)
#
#   And sets:
#
#     HAVE_BOOST
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2009 Peter Adolphs
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 45

# example boost program (need to pass version)
m4_define([_AX_BOOST_BASE_PROGRAM],
          [AC_LANG_PROGRAM([[
#include <boost/version.hpp>
]],[[
(void) ((void)sizeof(char[1 - 2*!!((BOOST_VERSION) < ($1))]));
]])])

AC_DEFUN([AX_BOOST_BASE],
[
AC_ARG_WITH([boost],
  [AS_HELP_STRING([--with-boost@<:@=ARG@:>@],
    [use Boost library from a standard location (ARG=yes),
     from the specified location (ARG=<path>),
     or disable it (ARG=no)
     @<:@ARG=yes@:>@ ])],
    [
     AS_CASE([$withval],
       [no],[want_boost="no";_AX_BOOST_BASE_boost_path=""],
       [yes],[want_boost="yes";_AX_BOOST_BASE_boost_path=""],
       [want_boost="yes";_AX_BOOST_BASE_boost_path="$withval"])
    ],
    [want_boost="yes"])


AC_ARG_WITH([boost-libdir],
  [AS_HELP_STRING([--with-boost-libdir=LIB_DIR],
    [Force given directory for boost libraries.
     Note that this will override library path detection,
     so use this parameter only if default library detection fails
     and you know exactly where your boost libraries are located.])],
  [
   AS_IF([test -d "$withval"],
         [_AX_BOOST_BASE_boost_lib_path="$withval"],
    [AC_MSG_ERROR([--with-boost-libdir expected directory name])])
  ],
  [_AX_BOOST_BASE_boost_lib_path=""])

BOOST_LDFLAGS=""
BOOST_CPPFLAGS=""
AS_IF([test "x$want_boost" = "xyes"],
      [_AX_BOOST_BASE_RUNDETECT([$1],[$2],[$3])])
AC_SUBST(BOOST_CPPFLAGS)
AC_SUBST(BOOST_LDFLAGS)
])


# convert a version string in $2 to numeric and affect to polymorphic var $1
AC_DEFUN([_AX_BOOST_BASE_TONUMERICVERSION],[
  AS_IF([test "x$2" = "x"],[_AX_BOOST_BASE_TONUMERICVERSION_req="1.20.0"],[_AX_BOOST_BASE_TONUMERICVERSION_req="$2"])
  _AX_BOOST_BASE_TONUMERICVERSION_req_shorten=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '\([[0-9]]*\.[[0-9]]*\)'`
  _AX_BOOST_BASE_TONUMERICVERSION_req_major=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '\([[0-9]]*\)'`
  AS_IF([test "x$_AX_BOOST_BASE_TONUMERICVERSION_req_major" = "x"],
        [AC_MSG_ERROR([You should at least specify libboost major version])])
  _AX_BOOST_BASE_TONUMERICVERSION_req_minor=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '[[0-9]]*\.\([[0-9]]*\)'`
  AS_IF([test "x$_AX_BOOST_BASE_TONUMERICVERSION_req_minor" = "x"],
        [_AX_BOOST_BASE_TONUMERICVERSION_req_minor="0"])
  _AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
  AS_IF([test "X$_AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor" = "X"],
        [_AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor="0"])
  _AX_BOOST_BASE_TONUMERICVERSION_RET=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req_major \* 100000 \+  $_AX_BOOST_BASE_TONUMERICVERSION_req_minor \* 100 \+ $_AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor`
  AS_VAR_SET($1,$_AX_BOOST_BASE_TONUMERICVERSION_RET)
])

dnl Run the detection of boost should be run only if $want_boost
AC_DEFUN([_AX_BOOST_BASE_RUNDETECT],[
    _AX_BOOST_BASE_TONUMERICVERSION(WANT_BOOST_VERSION,[$1])
    succeeded=no


    AC_REQUIRE([AC_CANONICAL_HOST])
    dnl On 64-bit systems check for system libraries in both lib64 and lib.
    dnl The former is specified by FHS, but e.g. Debian does not adhere to
    dnl this (as it rises problems for generic multi-arch support).
    dnl The last entry in the list is chosen by default when no libraries
    dnl are found, e.g. when only header-only libraries are installed!
    AS_CASE([${host_cpu}],
      [x86_64],[libsubdirs="lib64 libx32 lib lib64"],
      [ppc64|powerpc64|s390x|sparc64|aarch64|ppc64le|powerpc64le|riscv64],[libsubdirs="lib64 lib lib64"],
      [libsubdirs="lib"]
    )

    dnl allow for real multi-arch paths e.g. /usr/lib/x86_64-linux-gnu. Give
    dnl them priority over the other paths since, if libs are found there, they
    dnl are almost assuredly the ones desired.
    AS_CASE([${host_cpu}],
      [i?86],[multiarch_libsubdir="lib/i386-${host_os}"],
      [multiarch_libsubdir="lib/${host_cpu}-${host_os}"]
    )

    dnl first we check the system location for boost libraries
    dnl this location ist chosen if boost libraries are installed with the --layout=system option
    dnl or if you install boost with RPM
    AS_IF([test "x$_AX_BOOST_BASE_boost_path" != "x"],[
        AC_MSG_CHECKING([for boostlib >= $1 ($WANT_BOOST_VERSION) includes in "$_AX_BOOST_BASE_boost_path/include"])
         AS_IF([test -d "$_AX_BOOST_BASE_boost_path/include" && test -r "$_AX_BOOST_BASE_boost_path/include"],[
           AC_MSG_RESULT([yes])
           BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path/include"
           for _AX_BOOST_BASE_boost_path_tmp in $multiarch_libsubdir $libsubdirs; do
                AC_MSG_CHECKING([for boostlib >= $1 ($WANT_BOOST_VERSION) lib path in "$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp"])
                AS_IF([test -d "$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp" && test -r "$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp" ],[
                        AC_MSG_RESULT([yes])
                        BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp";
                        break;
                ],
      [AC_MSG_RESULT([no])])
           done],[
      AC_MSG_RESULT([no])])
    ],[
        if test X"$cross_compiling" = Xyes; then
            search_libsubdirs=$multiarch_libsubdir
        else
            search_libsubdirs="$multiarch_libsubdir $libsubdirs"
        fi
        for _AX_BOOST_BASE_boost_path_tmp in /usr /usr/local /opt /opt/local ; do
            if test -d "$_AX_BOOST_BASE_boost_path_tmp/include/boost" && test -r "$_AX_BOOST_BASE_boost_path_tmp/include/boost" ; then
                for libsubdir in $search_libsubdirs ; do
                    if ls "$_AX_BOOST_BASE_boost_path_tmp/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                done
                BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_path_tmp/$libsubdir"
                BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path_tmp/include"
                break;
            fi
        done
    ])

    dnl overwrite ld flags if we have required special directory with
    dnl --with-boost-libdir parameter
    AS_IF([test "x$_AX_BOOST_BASE_boost_lib_path" != "x"],
          [BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_lib_path"])

    AC_MSG_CHECKING([for boostlib >= $1 ($WANT_BOOST_VERSION)])
    CPPFLAGS_SAVED="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
    export CPPFLAGS

    LDFLAGS_SAVED="$LDFLAGS"
    LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
    export LDFLAGS

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)
        AC_COMPILE_IFELSE([_AX_BOOST_BASE_PROGRAM($WANT_BOOST_VERSION)],[
        AC_MSG_RESULT(yes)
    succeeded=yes
    found_system=yes
        ],[
        ])
    AC_LANG_POP([C++])



    dnl if we found no boost with system layout we search for boost libraries
    dnl built and installed without the --layout=system option or for a staged(not installed) version
    if test "x$succeeded" != "xyes" ; then
        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
        BOOST_CPPFLAGS=
        if test -z "$_AX_BOOST_BASE_boost_lib_path" ; then
            BOOST_LDFLAGS=
        fi
        _version=0
        if test -n "$_AX_BOOST_BASE_boost_path" ; then
            if test -d "$_AX_BOOST_BASE_boost_path" && test -r "$_AX_BOOST_BASE_boost_path"; then
                for i in `ls -d $_AX_BOOST_BASE_boost_path/include/boost-* 2>/dev/null`; do
                    _version_tmp=`echo $i | sed "s#$_AX_BOOST_BASE_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
                    V_CHECK=`expr $_version_tmp \> $_version`
                    if test "x$V_CHECK" = "x1" ; then
                        _version=$_version_tmp
                    fi
                    VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
                    BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path/include/boost-$VERSION_UNDERSCORE"
                done
                dnl if nothing found search for layout used in Windows distributions
                if test -z "$BOOST_CPPFLAGS"; then
                    if test -d "$_AX_BOOST_BASE_boost_path/boost" && test -r "$_AX_BOOST_BASE_boost_path/boost"; then
                        BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path"
                    fi
                fi
                dnl if we found something and BOOST_LDFLAGS was unset before
                dnl (because "$_AX_BOOST_BASE_boost_lib_path" = ""), set it here.
                if test -n "$BOOST_CPPFLAGS" && test -z "$BOOST_LDFLAGS"; then
                    for libsubdir in $libsubdirs ; do
                        if ls "$_AX_BOOST_BASE_boost_path/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                    done
                    BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_path/$libsubdir"
                fi
            fi
        else
            if test "x$cross_compiling" != "xyes" ; then
                for _AX_BOOST_BASE_boost_path in /usr /usr/local /opt /opt/local ; do
                    if test -d "$_AX_BOOST_BASE_boost_path" && test -r "$_AX_BOOST_BASE_boost_path" ; then
                        for i in `ls -d $_AX_BOOST_BASE_boost_path/include/boost-* 2>/dev/null`; do
                            _version_tmp=`echo $i | sed "s#$_AX_BOOST_BASE_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
                            V_CHECK=`expr $_version_tmp \> $_version`
                            if test "x$V_CHECK" = "x1" ; then
                                _version=$_version_tmp
                                best_path=$_AX_BOOST_BASE_boost_path
                            fi
                        done
                    fi
                done

                VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
                BOOST_CPPFLAGS="-I$best_path/include/boost-$VERSION_UNDERSCORE"
                if test -z "$_AX_BOOST_BASE_boost_lib_path" ; then
                    for libsubdir in $libsubdirs ; do
                        if ls "$best_path/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                    done
                    BOOST_LDFLAGS="-L$best_path/$libsubdir"
                fi
            fi

            if test -n "$BOOST_ROOT" ; then
                for libsubdir in $libsubdirs ; do
                    if ls "$BOOST_ROOT/stage/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                done
                if test -d "$BOOST_ROOT" && test -r "$BOOST_ROOT" && test -d "$BOOST_ROOT/stage/$libsubdir" && test -r "$BOOST_ROOT/stage/$libsubdir"; then
                    version_dir=`expr //$BOOST_ROOT : '.*/\(.*\)'`
                    stage_version=`echo $version_dir | sed 's/boost_//' | sed 's/_/./g'`
                        stage_version_shorten=`expr $stage_version : '\([[0-9]]*\.[[0-9]]*\)'`
                    V_CHECK=`expr $stage_version_shorten \>\= $_version`
                    if test "x$V_CHECK" = "x1" && test -z "$_AX_BOOST_BASE_boost_lib_path" ; then
                        AC_MSG_NOTICE(We will use a staged boost library from $BOOST_ROOT)
                        BOOST_CPPFLAGS="-I$BOOST_ROOT"
                        BOOST_LDFLAGS="-L$BOOST_ROOT/stage/$libsubdir"
                    fi
                fi
            fi
        fi

        CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
        export CPPFLAGS
        LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
        export LDFLAGS

        AC_LANG_PUSH(C++)
            AC_COMPILE_IFELSE([_AX_BOOST_BASE_PROGRAM($WANT_BOOST_VERSION)],[
            AC_MSG_RESULT(yes)
        succeeded=yes
        found_system=yes
            ],[
            ])
        AC_LANG_POP([C++])
    fi

    if test "x$succeeded" != "xyes" ; then
        if test "x$_version" = "x0" ; then
            AC_MSG_NOTICE([[We could not detect the boost libraries (version $1 or higher). If you have a staged boost library (still not installed) please specify \$BOOST_ROOT in your environment and do not give a PATH to --with-boost option.  If you are sure you have boost installed, then check your version number looking in <boost/version.hpp>. See http://randspringer.de/boost for more documentation.]])
        else
            AC_MSG_NOTICE([Your boost libraries seems to old (version $_version).])
        fi
        # execute ACTION-IF-NOT-FOUND (if present):
        ifelse([$3], , :, [$3])
    else
        AC_DEFINE(HAVE_BOOST,,[define if the Boost library is available])
        # execute ACTION-IF-FOUND (if present):
        ifelse([$2], , :, [$2])
    fi

    CPPFLAGS="$CPPFLAGS_SAVED"
    LDFLAGS="$LDFLAGS_SAVED"

])

# ===========================================================================
#    https://www.gnu.org/software/autoconf-archive/ax_boost_iostreams.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_IOSTREAMS
#
# DESCRIPTION
#
#   Test for IOStreams library from the Boost C++ libraries. The macro
#   requires a preceding call to AX_BOOST_BASE. Further documentation is
#   available at <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_IOSTREAMS_LIB)
#
#   And sets:
#
#     HAVE_BOOST_IOSTREAMS
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 21

AC_DEFUN([AX_BOOST_IOSTREAMS],
[
	AC_ARG_WITH([boost-iostreams],
	AS_HELP_STRING([--with-boost-iostreams@<:@=special-lib@:>@],
                   [use the IOStreams library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-iostreams=boost_iostreams-gcc-mt-d-1_33_1 ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_iostreams_lib=""
        else
		    want_boost="yes"
		ax_boost_user_iostreams_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::IOStreams library is available,
					   ax_cv_boost_iostreams,
        [AC_LANG_PUSH([C++])
		 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/iostreams/filtering_stream.hpp>
											 @%:@include <boost/range/iterator_range.hpp>
											]],
                                  [[std::string  input = "Hello World!";
								 namespace io = boost::iostreams;
									 io::filtering_istream  in(boost::make_iterator_range(input));
									 return 0;
                                   ]])],
                             ax_cv_boost_iostreams=yes, ax_cv_boost_iostreams=no)
         AC_LANG_POP([C++])
		])
		if test "x$ax_cv_boost_iostreams" = "xyes"; then
			AC_DEFINE(HAVE_BOOST_IOSTREAMS,,[define if the Boost::IOStreams library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
            if test "x$ax_boost_user_iostreams_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_iostreams*.so* $BOOSTLIBDIR/libboost_iostream*.dylib* $BOOSTLIBDIR/libboost_iostreams*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_iostreams.*\)\.so.*$;\1;' -e 's;^lib\(boost_iostream.*\)\.dylib.*$;\1;' -e 's;^lib\(boost_iostreams.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_IOSTREAMS_LIB="-l$ax_lib"; AC_SUBST(BOOST_IOSTREAMS_LIB) link_iostreams="yes"; break],
                                 [link_iostreams="no"])
				done
                if test "x$link_iostreams" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_iostreams*.dll* $BOOSTLIBDIR/boost_iostreams*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_iostreams.*\)\.dll.*$;\1;' -e 's;^\(boost_iostreams.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_IOSTREAMS_LIB="-l$ax_lib"; AC_SUBST(BOOST_IOSTREAMS_LIB) link_iostreams="yes"; break],
                                 [link_iostreams="no"])
				done
                fi

            else
               for ax_lib in $ax_boost_user_iostreams_lib boost_iostreams-$ax_boost_user_iostreams_lib; do
				      AC_CHECK_LIB($ax_lib, main,
                                   [BOOST_IOSTREAMS_LIB="-l$ax_lib"; AC_SUBST(BOOST_IOSTREAMS_LIB) link_iostreams="yes"; break],
                                   [link_iostreams="no"])
                  done

            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the library!)
            fi
			if test "x$link_iostreams" != "xyes"; then
				AC_MSG_ERROR(Could not link against $ax_lib !)
			fi
		fi

		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])

# =============================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_boost_program_options.html
# =============================================================================
#
# SYNOPSIS
#
#   AX_BOOST_PROGRAM_OPTIONS
#
# DESCRIPTION
#
#   Test for program options library from the Boost C++ libraries. The macro
#   requires a preceding call to AX_BOOST_BASE. Further documentation is
#   available at <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB)
#
#   And sets:
#
#     HAVE_BOOST_PROGRAM_OPTIONS
#
# LICENSE
#
#   Copyright (c) 2009 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 25

AC_DEFUN([AX_BOOST_PROGRAM_OPTIONS],
[
	AC_ARG_WITH([boost-program-options],
		AS_HELP_STRING([--with-boost-program-options@<:@=special-lib@:>@],
                       [use the program options library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-program-options=boost_program_options-gcc-mt-1_33_1 ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_program_options_lib=""
        else
		    want_boost="yes"
		ax_boost_user_program_options_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
	    export want_boost
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS
		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS
		AC_CACHE_CHECK([whether the Boost::Program_Options library is available],
					   ax_cv_boost_program_options,
					   [AC_LANG_PUSH(C++)
				AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/program_options/errors.hpp>
                                                          ]],
                                  [[boost::program_options::error err("Error message");
                                   return 0;]])],
                           ax_cv_boost_program_options=yes, ax_cv_boost_program_options=no)
					AC_LANG_POP([C++])
		])
		if test "$ax_cv_boost_program_options" = yes; then
				AC_DEFINE(HAVE_BOOST_PROGRAM_OPTIONS,,[define if the Boost::PROGRAM_OPTIONS library is available])
                  BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
                if test "x$ax_boost_user_program_options_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_program_options*.so* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_program_options.*\)\.so.*$;\1;'` `ls $BOOSTLIBDIR/libboost_program_options*.dylib* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_program_options.*\)\.dylib.*$;\1;'` `ls $BOOSTLIBDIR/libboost_program_options*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_program_options.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_PROGRAM_OPTIONS_LIB="-l$ax_lib"; AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB) link_program_options="yes"; break],
                                 [link_program_options="no"])
				done
                if test "x$link_program_options" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_program_options*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_program_options.*\)\.dll.*$;\1;'` `ls $BOOSTLIBDIR/boost_program_options*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_program_options.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_PROGRAM_OPTIONS_LIB="-l$ax_lib"; AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB) link_program_options="yes"; break],
                                 [link_program_options="no"])
				done
                fi
                else
                  for ax_lib in $ax_boost_user_program_options_lib boost_program_options-$ax_boost_user_program_options_lib; do
				      AC_CHECK_LIB($ax_lib, main,
                                   [BOOST_PROGRAM_OPTIONS_LIB="-l$ax_lib"; AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB) link_program_options="yes"; break],
                                   [link_program_options="no"])
                  done
                fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the library!)
            fi
				if test "x$link_program_options" != "xyes"; then
					AC_MSG_ERROR([Could not link against [$ax_lib] !])
				fi
		fi
		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])

# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_boost_regex.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_REGEX
#
# DESCRIPTION
#
#   Test for Regex library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_REGEX_LIB)
#
#   And sets:
#
#     HAVE_BOOST_REGEX
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Michael Tindal
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 23

AC_DEFUN([AX_BOOST_REGEX],
[
	AC_ARG_WITH([boost-regex],
	AS_HELP_STRING([--with-boost-regex@<:@=special-lib@:>@],
                   [use the Regex library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-regex=boost_regex-gcc-mt-d-1_33_1 ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_regex_lib=""
        else
		    want_boost="yes"
		ax_boost_user_regex_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::Regex library is available,
					   ax_cv_boost_regex,
        [AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/regex.hpp>
												]],
                                   [[boost::regex r(); return 0;]])],
                   ax_cv_boost_regex=yes, ax_cv_boost_regex=no)
         AC_LANG_POP([C++])
		])
		if test "x$ax_cv_boost_regex" = "xyes"; then
			AC_DEFINE(HAVE_BOOST_REGEX,,[define if the Boost::Regex library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
            if test "x$ax_boost_user_regex_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_regex*.so* $BOOSTLIBDIR/libboost_regex*.dylib* $BOOSTLIBDIR/libboost_regex*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_regex.*\)\.so.*$;\1;' -e 's;^lib\(boost_regex.*\)\.dylib.*;\1;' -e 's;^lib\(boost_regex.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_REGEX_LIB="-l$ax_lib"; AC_SUBST(BOOST_REGEX_LIB) link_regex="yes"; break],
                                 [link_regex="no"])
				done
                if test "x$link_regex" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_regex*.dll* $BOOSTLIBDIR/boost_regex*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_regex.*\)\.dll.*$;\1;' -e 's;^\(boost_regex.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_REGEX_LIB="-l$ax_lib"; AC_SUBST(BOOST_REGEX_LIB) link_regex="yes"; break],
                                 [link_regex="no"])
				done
                fi

            else
               for ax_lib in $ax_boost_user_regex_lib boost_regex-$ax_boost_user_regex_lib; do
				      AC_CHECK_LIB($ax_lib, main,
                                   [BOOST_REGEX_LIB="-l$ax_lib"; AC_SUBST(BOOST_REGEX_LIB) link_regex="yes"; break],
                                   [link_regex="no"])
               done
            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the Boost::Regex library!)
            fi
			if test "x$link_regex" != "xyes"; then
				AC_MSG_ERROR(Could not link against $ax_lib !)
			fi
		fi

		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])

# ===========================================================================
#     https://www.gnu.org/software/autoconf-archive/ax_boost_thread.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_THREAD
#
# DESCRIPTION
#
#   Test for Thread library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_THREAD_LIB)
#
#   And sets:
#
#     HAVE_BOOST_THREAD
#
# LICENSE
#
#   Copyright (c) 2009 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2009 Michael Tindal
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 31

AC_DEFUN([AX_BOOST_THREAD],
[
    AC_ARG_WITH([boost-thread],
    AS_HELP_STRING([--with-boost-thread@<:@=special-lib@:>@],
                   [use the Thread library from boost -
                    it is possible to specify a certain library for the linker
                    e.g. --with-boost-thread=boost_thread-gcc-mt ]),
        [
        if test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_thread_lib=""
        else
            want_boost="yes"
            ax_boost_user_thread_lib="$withval"
        fi
        ],
        [want_boost="yes"]
    )

    if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
        AC_REQUIRE([AC_CANONICAL_BUILD])
        CPPFLAGS_SAVED="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
        export CPPFLAGS

        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
        export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::Thread library is available,
                       ax_cv_boost_thread,
        [AC_LANG_PUSH([C++])
             CXXFLAGS_SAVE=$CXXFLAGS

             if test "x$host_os" = "xsolaris" ; then
                 CXXFLAGS="-pthreads $CXXFLAGS"
             elif test "x$host_os" = "xmingw32" ; then
                 CXXFLAGS="-mthreads $CXXFLAGS"
             else
                CXXFLAGS="-pthread $CXXFLAGS"
             fi
             AC_COMPILE_IFELSE([
                 AC_LANG_PROGRAM(
                     [[@%:@include <boost/thread/thread.hpp>]],
                     [[boost::thread_group thrds;
                       return 0;]])],
                 ax_cv_boost_thread=yes, ax_cv_boost_thread=no)
             CXXFLAGS=$CXXFLAGS_SAVE
             AC_LANG_POP([C++])
        ])
        if test "x$ax_cv_boost_thread" = "xyes"; then
           if test "x$host_os" = "xsolaris" ; then
              BOOST_CPPFLAGS="-pthreads $BOOST_CPPFLAGS"
           elif test "x$host_os" = "xmingw32" ; then
              BOOST_CPPFLAGS="-mthreads $BOOST_CPPFLAGS"
           else
              BOOST_CPPFLAGS="-pthread $BOOST_CPPFLAGS"
           fi

            AC_SUBST(BOOST_CPPFLAGS)

            AC_DEFINE(HAVE_BOOST_THREAD,,
                      [define if the Boost::Thread library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`

            LDFLAGS_SAVE=$LDFLAGS
                        case "x$host_os" in
                          *bsd* )
                               LDFLAGS="-pthread $LDFLAGS"
                          break;
                          ;;
                        esac
            if test "x$ax_boost_user_thread_lib" = "x"; then
                for libextension in `ls -r $BOOSTLIBDIR/libboost_thread* 2>/dev/null | sed 's,.*/lib,,' | sed 's,\..*,,'`; do
                     ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                                 [link_thread="yes"; break],
                                 [link_thread="no"])
                done
                if test "x$link_thread" != "xyes"; then
                for libextension in `ls -r $BOOSTLIBDIR/boost_thread* 2>/dev/null | sed 's,.*/,,' | sed 's,\..*,,'`; do
                     ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                                 [link_thread="yes"; break],
                                 [link_thread="no"])
                done
                fi

            else
               for ax_lib in $ax_boost_user_thread_lib boost_thread-$ax_boost_user_thread_lib; do
                      AC_CHECK_LIB($ax_lib, exit,
                                   [link_thread="yes"; break],
                                   [link_thread="no"])
                  done

            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the library!)
            fi
            if test "x$link_thread" = "xno"; then
                AC_MSG_ERROR(Could not link against $ax_lib !)
            else
                BOOST_THREAD_LIB="-l$ax_lib"
                case "x$host_os" in
                    *bsd* )
                        BOOST_LDFLAGS="-pthread $BOOST_LDFLAGS"
                        break;
                        ;;
                    xsolaris )
                        BOOST_THREAD_LIB="$BOOST_THREAD_LIB -lpthread"
                        break;
                        ;;
                    xmingw32 )
                        break;
                        ;;
                    * )
                        BOOST_THREAD_LIB="$BOOST_THREAD_LIB -lpthread"
                        break;
                        ;;
                esac
                AC_SUBST(BOOST_THREAD_LIB)
            fi
        fi

        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
    fi
])

# ===========================================================================
#     https://www.gnu.org/software/autoconf-archive/ax_check_library.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_LIBRARY(VARIABLE-PREFIX, HEADER-FILE, LIBRARY-FILE,
#                    [ACTION-IF-FOUND], [ACTION-IF-NOT_FOUND])
#
# DESCRIPTION
#
#   Provides a generic test for a given library, similar in concept to the
#   PKG_CHECK_MODULES macro used by pkg-config.
#
#   Most simplest libraries can be checked against simply through the
#   presence of a header file and a library to link to. This macro allows to
#   wrap around the test so that it doesn't have to be recreated each time.
#
#   Rather than define --with-$LIBRARY arguments, it uses variables in the
#   same way that PKG_CHECK_MODULES does. It doesn't, though, use the same
#   names, since you shouldn't provide a value for LIBS or CFLAGS but rather
#   for LDFLAGS and CPPFLAGS, to tell the linker and compiler where to find
#   libraries and headers respectively.
#
#   If the library is find, HAVE_PREFIX is defined, and in all cases
#   PREFIX_LDFLAGS and PREFIX_CPPFLAGS are substituted.
#
#   Example:
#
#     AX_CHECK_LIBRARY([LIBEVENT], [event.h], [event], [],
#                      [AC_MSG_ERROR([Unable to find libevent])])
#
# LICENSE
#
#   Copyright (c) 2010 Diego Elio Petteno` <flameeyes@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 5

AC_DEFUN([AX_CHECK_LIBRARY], [
  AC_ARG_VAR($1[_CPPFLAGS], [C preprocessor flags for ]$1[ headers])
  AC_ARG_VAR($1[_LDFLAGS], [linker flags for ]$1[ libraries])

  AC_CACHE_VAL(AS_TR_SH([ax_cv_have_]$1),
    [save_CPPFLAGS="$CPPFLAGS"
     save_LDFLAGS="$LDFLAGS"
     save_LIBS="$LIBS"

     AS_IF([test "x$]$1[_CPPFLAGS" != "x"],
       [CPPFLAGS="$CPPFLAGS $]$1[_CPPFLAGS"])

     AS_IF([test "x$]$1[_LDFLAGS" != "x"],
       [LDFLAGS="$LDFLAGS $]$1[_LDFLAGS"])

     AC_CHECK_HEADER($2, [
       AC_CHECK_LIB($3, [main],
         [AS_TR_SH([ax_cv_have_]$1)=yes],
         [AS_TR_SH([ax_cv_have_]$1)=no])
     ], [AS_TR_SH([ax_cv_have_]$1)=no])

     CPPFLAGS="$save_CPPFLAGS"
     LDFLAGS="$save_LDFLAGS"
     LIBS="$save_LIBS"
    ])

  AS_IF([test "$]AS_TR_SH([ax_cv_have_]$1)[" = "yes"],
    AC_DEFINE([HAVE_]$1, [1], [Define to 1 if ]$1[ is found])
    [$4],
    [$5])
])

# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX(VERSION, [ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the specified
#   version of the C++ standard.  If necessary, add switches to CXX and
#   CXXCPP to enable support.  VERSION may be '11' (for the C++11 standard)
#   or '14' (for the C++14 standard).
#
#   The second argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for an extended mode.
#
#   The third argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline support for the specified C++ standard is
#   required and that the macro should error out if no mode with that
#   support is found.  If specified 'optional', then configuration proceeds
#   regardless, after defining HAVE_CXX${VERSION} if and only if a
#   supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014, 2015 Google Inc.; contributed by Alexey Sokolov <sokolov@google.com>
#   Copyright (c) 2015 Paul Norman <penorman@mac.com>
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#   Copyright (c) 2016, 2018 Krzesimir Nowak <qdlacz@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 10

dnl  This macro is based on the code from the AX_CXX_COMPILE_STDCXX_11 macro
dnl  (serial version number 13).

AC_DEFUN([AX_CXX_COMPILE_STDCXX], [dnl
  m4_if([$1], [11], [ax_cxx_compile_alternatives="11 0x"],
        [$1], [14], [ax_cxx_compile_alternatives="14 1y"],
        [$1], [17], [ax_cxx_compile_alternatives="17 1z"],
        [m4_fatal([invalid first argument `$1' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$2], [], [],
        [$2], [ext], [],
        [$2], [noext], [],
        [m4_fatal([invalid second argument `$2' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$3], [], [ax_cxx_compile_cxx$1_required=true],
        [$3], [mandatory], [ax_cxx_compile_cxx$1_required=true],
        [$3], [optional], [ax_cxx_compile_cxx$1_required=false],
        [m4_fatal([invalid third argument `$3' to AX_CXX_COMPILE_STDCXX])])
  AC_LANG_PUSH([C++])dnl
  ac_success=no

  m4_if([$2], [noext], [], [dnl
  if test x$ac_success = xno; then
    for alternative in ${ax_cxx_compile_alternatives}; do
      switch="-std=gnu++${alternative}"
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                     $cachevar,
        [ac_save_CXX="$CXX"
         CXX="$CXX $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXX="$ac_save_CXX"])
      if eval test x\$$cachevar = xyes; then
        CXX="$CXX $switch"
        if test -n "$CXXCPP" ; then
          CXXCPP="$CXXCPP $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])

  m4_if([$2], [ext], [], [dnl
  if test x$ac_success = xno; then
    dnl HP's aCC needs +std=c++11 according to:
    dnl http://h21007.www2.hp.com/portal/download/files/unprot/aCxx/PDF_Release_Notes/769149-001.pdf
    dnl Cray's crayCC needs "-h std=c++11"
    for alternative in ${ax_cxx_compile_alternatives}; do
      for switch in -std=c++${alternative} +std=c++${alternative} "-h std=c++${alternative}"; do
        cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
        AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                       $cachevar,
          [ac_save_CXX="$CXX"
           CXX="$CXX $switch"
           AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
            [eval $cachevar=yes],
            [eval $cachevar=no])
           CXX="$ac_save_CXX"])
        if eval test x\$$cachevar = xyes; then
          CXX="$CXX $switch"
          if test -n "$CXXCPP" ; then
            CXXCPP="$CXXCPP $switch"
          fi
          ac_success=yes
          break
        fi
      done
      if test x$ac_success = xyes; then
        break
      fi
    done
  fi])
  AC_LANG_POP([C++])
  if test x$ax_cxx_compile_cxx$1_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler with support for C++$1 language features is required.])
    fi
  fi
  if test x$ac_success = xno; then
    HAVE_CXX$1=0
    AC_MSG_NOTICE([No compiler with C++$1 support was found])
  else
    HAVE_CXX$1=1
    AC_DEFINE(HAVE_CXX$1,1,
              [define if the compiler supports basic C++$1 syntax])
  fi
  AC_SUBST(HAVE_CXX$1)
])


dnl  Test body for checking C++11 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_11],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
)


dnl  Test body for checking C++14 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_14],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
)

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_17],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_17
)

dnl  Tests for new features in C++11

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_11], [[

// If the compiler admits that it is not ready for C++11, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201103L

#error "This is not a C++11 compiler"

#else

namespace cxx11
{

  namespace test_static_assert
  {

    template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

  }

  namespace test_final_override
  {

    struct Base
    {
      virtual void f() {}
    };

    struct Derived : public Base
    {
      virtual void f() override {}
    };

  }

  namespace test_double_right_angle_brackets
  {

    template < typename T >
    struct check {};

    typedef check<void> single_type;
    typedef check<check<void>> double_type;
    typedef check<check<check<void>>> triple_type;
    typedef check<check<check<check<void>>>> quadruple_type;

  }

  namespace test_decltype
  {

    int
    f()
    {
      int a = 1;
      decltype(a) b = 2;
      return a + b;
    }

  }

  namespace test_type_deduction
  {

    template < typename T1, typename T2 >
    struct is_same
    {
      static const bool value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static const bool value = true;
    };

    template < typename T1, typename T2 >
    auto
    add(T1 a1, T2 a2) -> decltype(a1 + a2)
    {
      return a1 + a2;
    }

    int
    test(const int c, volatile int v)
    {
      static_assert(is_same<int, decltype(0)>::value == true, "");
      static_assert(is_same<int, decltype(c)>::value == false, "");
      static_assert(is_same<int, decltype(v)>::value == false, "");
      auto ac = c;
      auto av = v;
      auto sumi = ac + av + 'x';
      auto sumf = ac + av + 1.0;
      static_assert(is_same<int, decltype(ac)>::value == true, "");
      static_assert(is_same<int, decltype(av)>::value == true, "");
      static_assert(is_same<int, decltype(sumi)>::value == true, "");
      static_assert(is_same<int, decltype(sumf)>::value == false, "");
      static_assert(is_same<int, decltype(add(c, v))>::value == true, "");
      return (sumf > 0.0) ? sumi : add(c, v);
    }

  }

  namespace test_noexcept
  {

    int f() { return 0; }
    int g() noexcept { return 0; }

    static_assert(noexcept(f()) == false, "");
    static_assert(noexcept(g()) == true, "");

  }

  namespace test_constexpr
  {

    template < typename CharT >
    unsigned long constexpr
    strlen_c_r(const CharT *const s, const unsigned long acc) noexcept
    {
      return *s ? strlen_c_r(s + 1, acc + 1) : acc;
    }

    template < typename CharT >
    unsigned long constexpr
    strlen_c(const CharT *const s) noexcept
    {
      return strlen_c_r(s, 0UL);
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("1") == 1UL, "");
    static_assert(strlen_c("example") == 7UL, "");
    static_assert(strlen_c("another\0example") == 7UL, "");

  }

  namespace test_rvalue_references
  {

    template < int N >
    struct answer
    {
      static constexpr int value = N;
    };

    answer<1> f(int&)       { return answer<1>(); }
    answer<2> f(const int&) { return answer<2>(); }
    answer<3> f(int&&)      { return answer<3>(); }

    void
    test()
    {
      int i = 0;
      const int c = 0;
      static_assert(decltype(f(i))::value == 1, "");
      static_assert(decltype(f(c))::value == 2, "");
      static_assert(decltype(f(0))::value == 3, "");
    }

  }

  namespace test_uniform_initialization
  {

    struct test
    {
      static const int zero {};
      static const int one {1};
    };

    static_assert(test::zero == 0, "");
    static_assert(test::one == 1, "");

  }

  namespace test_lambdas
  {

    void
    test1()
    {
      auto lambda1 = [](){};
      auto lambda2 = lambda1;
      lambda1();
      lambda2();
    }

    int
    test2()
    {
      auto a = [](int i, int j){ return i + j; }(1, 2);
      auto b = []() -> int { return '0'; }();
      auto c = [=](){ return a + b; }();
      auto d = [&](){ return c; }();
      auto e = [a, &b](int x) mutable {
        const auto identity = [](int y){ return y; };
        for (auto i = 0; i < a; ++i)
          a += b--;
        return x + identity(a + b);
      }(0);
      return a + b + c + d + e;
    }

    int
    test3()
    {
      const auto nullary = [](){ return 0; };
      const auto unary = [](int x){ return x; };
      using nullary_t = decltype(nullary);
      using unary_t = decltype(unary);
      const auto higher1st = [](nullary_t f){ return f(); };
      const auto higher2nd = [unary](nullary_t f1){
        return [unary, f1](unary_t f2){ return f2(unary(f1())); };
      };
      return higher1st(nullary) + higher2nd(nullary)(unary);
    }

  }

  namespace test_variadic_templates
  {

    template <int...>
    struct sum;

    template <int N0, int... N1toN>
    struct sum<N0, N1toN...>
    {
      static constexpr auto value = N0 + sum<N1toN...>::value;
    };

    template <>
    struct sum<>
    {
      static constexpr auto value = 0;
    };

    static_assert(sum<>::value == 0, "");
    static_assert(sum<1>::value == 1, "");
    static_assert(sum<23>::value == 23, "");
    static_assert(sum<1, 2>::value == 3, "");
    static_assert(sum<5, 5, 11>::value == 21, "");
    static_assert(sum<2, 3, 5, 7, 11, 13>::value == 41, "");

  }

  // http://stackoverflow.com/questions/13728184/template-aliases-and-sfinae
  // Clang 3.1 fails with headers of libstd++ 4.8.3 when using std::function
  // because of this.
  namespace test_template_alias_sfinae
  {

    struct foo {};

    template<typename T>
    using member = typename T::member_type;

    template<typename T>
    void func(...) {}

    template<typename T>
    void func(member<T>*) {}

    void test();

    void test() { func<foo>(0); }

  }

}  // namespace cxx11

#endif  // __cplusplus >= 201103L

]])


dnl  Tests for new features in C++14

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_14], [[

// If the compiler admits that it is not ready for C++14, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201402L

#error "This is not a C++14 compiler"

#else

namespace cxx14
{

  namespace test_polymorphic_lambdas
  {

    int
    test()
    {
      const auto lambda = [](auto&&... args){
        const auto istiny = [](auto x){
          return (sizeof(x) == 1UL) ? 1 : 0;
        };
        const int aretiny[] = { istiny(args)... };
        return aretiny[0];
      };
      return lambda(1, 1L, 1.0f, '1');
    }

  }

  namespace test_binary_literals
  {

    constexpr auto ivii = 0b0000000000101010;
    static_assert(ivii == 42, "wrong value");

  }

  namespace test_generalized_constexpr
  {

    template < typename CharT >
    constexpr unsigned long
    strlen_c(const CharT *const s) noexcept
    {
      auto length = 0UL;
      for (auto p = s; *p; ++p)
        ++length;
      return length;
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("x") == 1UL, "");
    static_assert(strlen_c("test") == 4UL, "");
    static_assert(strlen_c("another\0test") == 7UL, "");

  }

  namespace test_lambda_init_capture
  {

    int
    test()
    {
      auto x = 0;
      const auto lambda1 = [a = x](int b){ return a + b; };
      const auto lambda2 = [a = lambda1(x)](){ return a; };
      return lambda2();
    }

  }

  namespace test_digit_separators
  {

    constexpr auto ten_million = 100'000'000;
    static_assert(ten_million == 100000000, "");

  }

  namespace test_return_type_deduction
  {

    auto f(int& x) { return x; }
    decltype(auto) g(int& x) { return x; }

    template < typename T1, typename T2 >
    struct is_same
    {
      static constexpr auto value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static constexpr auto value = true;
    };

    int
    test()
    {
      auto x = 0;
      static_assert(is_same<int, decltype(f(x))>::value, "");
      static_assert(is_same<int&, decltype(g(x))>::value, "");
      return x;
    }

  }

}  // namespace cxx14

#endif  // __cplusplus >= 201402L

]])


dnl  Tests for new features in C++17

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_17], [[

// If the compiler admits that it is not ready for C++17, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201703L

#error "This is not a C++17 compiler"

#else

#include <initializer_list>
#include <utility>
#include <type_traits>

namespace cxx17
{

  namespace test_constexpr_lambdas
  {

    constexpr int foo = [](){return 42;}();

  }

  namespace test::nested_namespace::definitions
  {

  }

  namespace test_fold_expression
  {

    template<typename... Args>
    int multiply(Args... args)
    {
      return (args * ... * 1);
    }

    template<typename... Args>
    bool all(Args... args)
    {
      return (args && ...);
    }

  }

  namespace test_extended_static_assert
  {

    static_assert (true);

  }

  namespace test_auto_brace_init_list
  {

    auto foo = {5};
    auto bar {5};

    static_assert(std::is_same<std::initializer_list<int>, decltype(foo)>::value);
    static_assert(std::is_same<int, decltype(bar)>::value);
  }

  namespace test_typename_in_template_template_parameter
  {

    template<template<typename> typename X> struct D;

  }

  namespace test_fallthrough_nodiscard_maybe_unused_attributes
  {

    int f1()
    {
      return 42;
    }

    [[nodiscard]] int f2()
    {
      [[maybe_unused]] auto unused = f1();

      switch (f1())
      {
      case 17:
        f1();
        [[fallthrough]];
      case 42:
        f1();
      }
      return f1();
    }

  }

  namespace test_extended_aggregate_initialization
  {

    struct base1
    {
      int b1, b2 = 42;
    };

    struct base2
    {
      base2() {
        b3 = 42;
      }
      int b3;
    };

    struct derived : base1, base2
    {
        int d;
    };

    derived d1 {{1, 2}, {}, 4};  // full initialization
    derived d2 {{}, {}, 4};      // value-initialized bases

  }

  namespace test_general_range_based_for_loop
  {

    struct iter
    {
      int i;

      int& operator* ()
      {
        return i;
      }

      const int& operator* () const
      {
        return i;
      }

      iter& operator++()
      {
        ++i;
        return *this;
      }
    };

    struct sentinel
    {
      int i;
    };

    bool operator== (const iter& i, const sentinel& s)
    {
      return i.i == s.i;
    }

    bool operator!= (const iter& i, const sentinel& s)
    {
      return !(i == s);
    }

    struct range
    {
      iter begin() const
      {
        return {0};
      }

      sentinel end() const
      {
        return {5};
      }
    };

    void f()
    {
      range r {};

      for (auto i : r)
      {
        [[maybe_unused]] auto v = i;
      }
    }

  }

  namespace test_lambda_capture_asterisk_this_by_value
  {

    struct t
    {
      int i;
      int foo()
      {
        return [*this]()
        {
          return i;
        }();
      }
    };

  }

  namespace test_enum_class_construction
  {

    enum class byte : unsigned char
    {};

    byte foo {42};

  }

  namespace test_constexpr_if
  {

    template <bool cond>
    int f ()
    {
      if constexpr(cond)
      {
        return 13;
      }
      else
      {
        return 42;
      }
    }

  }

  namespace test_selection_statement_with_initializer
  {

    int f()
    {
      return 13;
    }

    int f2()
    {
      if (auto i = f(); i > 0)
      {
        return 3;
      }

      switch (auto i = f(); i + 4)
      {
      case 17:
        return 2;

      default:
        return 1;
      }
    }

  }

  namespace test_template_argument_deduction_for_class_templates
  {

    template <typename T1, typename T2>
    struct pair
    {
      pair (T1 p1, T2 p2)
        : m1 {p1},
          m2 {p2}
      {}

      T1 m1;
      T2 m2;
    };

    void f()
    {
      [[maybe_unused]] auto p = pair{13, 42u};
    }

  }

  namespace test_non_type_auto_template_parameters
  {

    template <auto n>
    struct B
    {};

    B<5> b1;
    B<'a'> b2;

  }

  namespace test_structured_bindings
  {

    int arr[2] = { 1, 2 };
    std::pair<int, int> pr = { 1, 2 };

    auto f1() -> int(&)[2]
    {
      return arr;
    }

    auto f2() -> std::pair<int, int>&
    {
      return pr;
    }

    struct S
    {
      int x1 : 2;
      volatile double y1;
    };

    S f3()
    {
      return {};
    }

    auto [ x1, y1 ] = f1();
    auto& [ xr1, yr1 ] = f1();
    auto [ x2, y2 ] = f2();
    auto& [ xr2, yr2 ] = f2();
    const auto [ x3, y3 ] = f3();

  }

  namespace test_exception_spec_type_system
  {

    struct Good {};
    struct Bad {};

    void g1() noexcept;
    void g2();

    template<typename T>
    Bad
    f(T*, T*);

    template<typename T1, typename T2>
    Good
    f(T1*, T2*);

    static_assert (std::is_same_v<Good, decltype(f(g1, g2))>);

  }

  namespace test_inline_variables
  {

    template<class T> void f(T)
    {}

    template<class T> inline T g(T)
    {
      return T{};
    }

    template<> inline void f<>(int)
    {}

    template<> int g<>(int)
    {
      return 5;
    }

  }

}  // namespace cxx17

#endif  // __cplusplus < 201703L

]])

# =============================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_17.html
# =============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_17([ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++17
#   standard; if necessary, add switches to CXX and CXXCPP to enable
#   support.
#
#   This macro is a convenience alias for calling the AX_CXX_COMPILE_STDCXX
#   macro with the version set to C++17.  The two optional arguments are
#   forwarded literally as the second and third argument respectively.
#   Please see the documentation for the AX_CXX_COMPILE_STDCXX macro for
#   more information.  If you want to use this macro, you also need to
#   download the ax_cxx_compile_stdcxx.m4 file.
#
# LICENSE
#
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#   Copyright (c) 2016 Krzesimir Nowak <qdlacz@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 2

AX_REQUIRE_DEFINED([AX_CXX_COMPILE_STDCXX])
AC_DEFUN([AX_CXX_COMPILE_STDCXX_17], [AX_CXX_COMPILE_STDCXX([17], [$1], [$2])])

# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_pthread.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro figures out how to build C programs using POSIX threads. It
#   sets the PTHREAD_LIBS output variable to the threads library and linker
#   flags, and the PTHREAD_CFLAGS output variable to any special C compiler
#   flags that are needed. (The user can also force certain compiler
#   flags/libs to be tested by setting these environment variables.)
#
#   Also sets PTHREAD_CC to any special C compiler that is needed for
#   multi-threaded programs (defaults to the value of CC otherwise). (This
#   is necessary on AIX to use the special cc_r compiler alias.)
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also to link with them as well. For example, you might link with
#   $PTHREAD_CC $CFLAGS $PTHREAD_CFLAGS $LDFLAGS ... $PTHREAD_LIBS $LIBS
#
#   If you are only building threaded programs, you may wish to use these
#   variables in your default LIBS, CFLAGS, and CC:
#
#     LIBS="$PTHREAD_LIBS $LIBS"
#     CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
#     CC="$PTHREAD_CC"
#
#   In addition, if the PTHREAD_CREATE_JOINABLE thread-attribute constant
#   has a nonstandard name, this macro defines PTHREAD_CREATE_JOINABLE to
#   that name (e.g. PTHREAD_CREATE_UNDETACHED on AIX).
#
#   Also HAVE_PTHREAD_PRIO_INHERIT is defined if pthread is found and the
#   PTHREAD_PRIO_INHERIT symbol is defined when compiling with
#   PTHREAD_CFLAGS.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a threads library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_PTHREAD.
#
#   Please let the authors know if this macro fails on any platform, or if
#   you have any other suggestions or comments. This macro was based on work
#   by SGJ on autoconf scripts for FFTW (http://www.fftw.org/) (with help
#   from M. Frigo), as well as ac_pthread and hb_pthread macros posted by
#   Alejandro Forero Cuervo to the autoconf macro repository. We are also
#   grateful for the helpful feedback of numerous users.
#
#   Updated for Autoconf 2.68 by Daniel Richard G.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2011 Daniel Richard G. <skunk@iSKUNK.ORG>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 24

AU_ALIAS([ACX_PTHREAD], [AX_PTHREAD])
AC_DEFUN([AX_PTHREAD], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_SED])
AC_LANG_PUSH([C])
ax_pthread_ok=no

# We used to check for pthread.h first, but this fails if pthread.h
# requires special compiler flags (e.g. on Tru64 or Sequent).
# It gets checked for in the link test anyway.

# First of all, check if the user has set any of the PTHREAD_LIBS,
# etcetera environment variables, and if threads linking works using
# them:
if test "x$PTHREAD_CFLAGS$PTHREAD_LIBS" != "x"; then
        ax_pthread_save_CC="$CC"
        ax_pthread_save_CFLAGS="$CFLAGS"
        ax_pthread_save_LIBS="$LIBS"
        AS_IF([test "x$PTHREAD_CC" != "x"], [CC="$PTHREAD_CC"])
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"
        AC_MSG_CHECKING([for pthread_join using $CC $PTHREAD_CFLAGS $PTHREAD_LIBS])
        AC_LINK_IFELSE([AC_LANG_CALL([], [pthread_join])], [ax_pthread_ok=yes])
        AC_MSG_RESULT([$ax_pthread_ok])
        if test "x$ax_pthread_ok" = "xno"; then
                PTHREAD_LIBS=""
                PTHREAD_CFLAGS=""
        fi
        CC="$ax_pthread_save_CC"
        CFLAGS="$ax_pthread_save_CFLAGS"
        LIBS="$ax_pthread_save_LIBS"
fi

# We must check for the threads library under a number of different
# names; the ordering is very important because some systems
# (e.g. DEC) have both -lpthread and -lpthreads, where one of the
# libraries is broken (non-POSIX).

# Create a list of thread flags to try.  Items starting with a "-" are
# C compiler flags, and other items are library names, except for "none"
# which indicates that we try without any flags at all, and "pthread-config"
# which is a program returning the flags for the Pth emulation library.

ax_pthread_flags="pthreads none -Kthread -pthread -pthreads -mthreads pthread --thread-safe -mt pthread-config"

# The ordering *is* (sometimes) important.  Some notes on the
# individual items follow:

# pthreads: AIX (must check this before -lpthread)
# none: in case threads are in libc; should be tried before -Kthread and
#       other compiler flags to prevent continual compiler warnings
# -Kthread: Sequent (threads in libc, but -Kthread needed for pthread.h)
# -pthread: Linux/gcc (kernel threads), BSD/gcc (userland threads), Tru64
#           (Note: HP C rejects this with "bad form for `-t' option")
# -pthreads: Solaris/gcc (Note: HP C also rejects)
# -mt: Sun Workshop C (may only link SunOS threads [-lthread], but it
#      doesn't hurt to check since this sometimes defines pthreads and
#      -D_REENTRANT too), HP C (must be checked before -lpthread, which
#      is present but should not be used directly; and before -mthreads,
#      because the compiler interprets this as "-mt" + "-hreads")
# -mthreads: Mingw32/gcc, Lynx/gcc
# pthread: Linux, etcetera
# --thread-safe: KAI C++
# pthread-config: use pthread-config program (for GNU Pth library)

case $host_os in

        freebsd*)

        # -kthread: FreeBSD kernel threads (preferred to -pthread since SMP-able)
        # lthread: LinuxThreads port on FreeBSD (also preferred to -pthread)

        ax_pthread_flags="-kthread lthread $ax_pthread_flags"
        ;;

        hpux*)

        # From the cc(1) man page: "[-mt] Sets various -D flags to enable
        # multi-threading and also sets -lpthread."

        ax_pthread_flags="-mt -pthread pthread $ax_pthread_flags"
        ;;

        openedition*)

        # IBM z/OS requires a feature-test macro to be defined in order to
        # enable POSIX threads at all, so give the user a hint if this is
        # not set. (We don't define these ourselves, as they can affect
        # other portions of the system API in unpredictable ways.)

        AC_EGREP_CPP([AX_PTHREAD_ZOS_MISSING],
            [
#            if !defined(_OPEN_THREADS) && !defined(_UNIX03_THREADS)
             AX_PTHREAD_ZOS_MISSING
#            endif
            ],
            [AC_MSG_WARN([IBM z/OS requires -D_OPEN_THREADS or -D_UNIX03_THREADS to enable pthreads support.])])
        ;;

        solaris*)

        # On Solaris (at least, for some versions), libc contains stubbed
        # (non-functional) versions of the pthreads routines, so link-based
        # tests will erroneously succeed. (N.B.: The stubs are missing
        # pthread_cleanup_push, or rather a function called by this macro,
        # so we could check for that, but who knows whether they'll stub
        # that too in a future libc.)  So we'll check first for the
        # standard Solaris way of linking pthreads (-mt -lpthread).

        ax_pthread_flags="-mt,pthread pthread $ax_pthread_flags"
        ;;
esac

# GCC generally uses -pthread, or -pthreads on some platforms (e.g. SPARC)

AS_IF([test "x$GCC" = "xyes"],
      [ax_pthread_flags="-pthread -pthreads $ax_pthread_flags"])

# The presence of a feature test macro requesting re-entrant function
# definitions is, on some systems, a strong hint that pthreads support is
# correctly enabled

case $host_os in
        darwin* | hpux* | linux* | osf* | solaris*)
        ax_pthread_check_macro="_REENTRANT"
        ;;

        aix*)
        ax_pthread_check_macro="_THREAD_SAFE"
        ;;

        *)
        ax_pthread_check_macro="--"
        ;;
esac
AS_IF([test "x$ax_pthread_check_macro" = "x--"],
      [ax_pthread_check_cond=0],
      [ax_pthread_check_cond="!defined($ax_pthread_check_macro)"])

# Are we compiling with Clang?

AC_CACHE_CHECK([whether $CC is Clang],
    [ax_cv_PTHREAD_CLANG],
    [ax_cv_PTHREAD_CLANG=no
     # Note that Autoconf sets GCC=yes for Clang as well as GCC
     if test "x$GCC" = "xyes"; then
        AC_EGREP_CPP([AX_PTHREAD_CC_IS_CLANG],
            [/* Note: Clang 2.7 lacks __clang_[a-z]+__ */
#            if defined(__clang__) && defined(__llvm__)
             AX_PTHREAD_CC_IS_CLANG
#            endif
            ],
            [ax_cv_PTHREAD_CLANG=yes])
     fi
    ])
ax_pthread_clang="$ax_cv_PTHREAD_CLANG"

ax_pthread_clang_warning=no

# Clang needs special handling, because older versions handle the -pthread
# option in a rather... idiosyncratic way

if test "x$ax_pthread_clang" = "xyes"; then

        # Clang takes -pthread; it has never supported any other flag

        # (Note 1: This will need to be revisited if a system that Clang
        # supports has POSIX threads in a separate library.  This tends not
        # to be the way of modern systems, but it's conceivable.)

        # (Note 2: On some systems, notably Darwin, -pthread is not needed
        # to get POSIX threads support; the API is always present and
        # active.  We could reasonably leave PTHREAD_CFLAGS empty.  But
        # -pthread does define _REENTRANT, and while the Darwin headers
        # ignore this macro, third-party headers might not.)

        PTHREAD_CFLAGS="-pthread"
        PTHREAD_LIBS=

        ax_pthread_ok=yes

        # However, older versions of Clang make a point of warning the user
        # that, in an invocation where only linking and no compilation is
        # taking place, the -pthread option has no effect ("argument unused
        # during compilation").  They expect -pthread to be passed in only
        # when source code is being compiled.
        #
        # Problem is, this is at odds with the way Automake and most other
        # C build frameworks function, which is that the same flags used in
        # compilation (CFLAGS) are also used in linking.  Many systems
        # supported by AX_PTHREAD require exactly this for POSIX threads
        # support, and in fact it is often not straightforward to specify a
        # flag that is used only in the compilation phase and not in
        # linking.  Such a scenario is extremely rare in practice.
        #
        # Even though use of the -pthread flag in linking would only print
        # a warning, this can be a nuisance for well-run software projects
        # that build with -Werror.  So if the active version of Clang has
        # this misfeature, we search for an option to squash it.

        AC_CACHE_CHECK([whether Clang needs flag to prevent "argument unused" warning when linking with -pthread],
            [ax_cv_PTHREAD_CLANG_NO_WARN_FLAG],
            [ax_cv_PTHREAD_CLANG_NO_WARN_FLAG=unknown
             # Create an alternate version of $ac_link that compiles and
             # links in two steps (.c -> .o, .o -> exe) instead of one
             # (.c -> exe), because the warning occurs only in the second
             # step
             ax_pthread_save_ac_link="$ac_link"
             ax_pthread_sed='s/conftest\.\$ac_ext/conftest.$ac_objext/g'
             ax_pthread_link_step=`$as_echo "$ac_link" | sed "$ax_pthread_sed"`
             ax_pthread_2step_ac_link="($ac_compile) && (echo ==== >&5) && ($ax_pthread_link_step)"
             ax_pthread_save_CFLAGS="$CFLAGS"
             for ax_pthread_try in '' -Qunused-arguments -Wno-unused-command-line-argument unknown; do
                AS_IF([test "x$ax_pthread_try" = "xunknown"], [break])
                CFLAGS="-Werror -Wunknown-warning-option $ax_pthread_try -pthread $ax_pthread_save_CFLAGS"
                ac_link="$ax_pthread_save_ac_link"
                AC_LINK_IFELSE([AC_LANG_SOURCE([[int main(void){return 0;}]])],
                    [ac_link="$ax_pthread_2step_ac_link"
                     AC_LINK_IFELSE([AC_LANG_SOURCE([[int main(void){return 0;}]])],
                         [break])
                    ])
             done
             ac_link="$ax_pthread_save_ac_link"
             CFLAGS="$ax_pthread_save_CFLAGS"
             AS_IF([test "x$ax_pthread_try" = "x"], [ax_pthread_try=no])
             ax_cv_PTHREAD_CLANG_NO_WARN_FLAG="$ax_pthread_try"
            ])

        case "$ax_cv_PTHREAD_CLANG_NO_WARN_FLAG" in
                no | unknown) ;;
                *) PTHREAD_CFLAGS="$ax_cv_PTHREAD_CLANG_NO_WARN_FLAG $PTHREAD_CFLAGS" ;;
        esac

fi # $ax_pthread_clang = yes

if test "x$ax_pthread_ok" = "xno"; then
for ax_pthread_try_flag in $ax_pthread_flags; do

        case $ax_pthread_try_flag in
                none)
                AC_MSG_CHECKING([whether pthreads work without any flags])
                ;;

                -mt,pthread)
                AC_MSG_CHECKING([whether pthreads work with -mt -lpthread])
                PTHREAD_CFLAGS="-mt"
                PTHREAD_LIBS="-lpthread"
                ;;

                -*)
                AC_MSG_CHECKING([whether pthreads work with $ax_pthread_try_flag])
                PTHREAD_CFLAGS="$ax_pthread_try_flag"
                ;;

                pthread-config)
                AC_CHECK_PROG([ax_pthread_config], [pthread-config], [yes], [no])
                AS_IF([test "x$ax_pthread_config" = "xno"], [continue])
                PTHREAD_CFLAGS="`pthread-config --cflags`"
                PTHREAD_LIBS="`pthread-config --ldflags` `pthread-config --libs`"
                ;;

                *)
                AC_MSG_CHECKING([for the pthreads library -l$ax_pthread_try_flag])
                PTHREAD_LIBS="-l$ax_pthread_try_flag"
                ;;
        esac

        ax_pthread_save_CFLAGS="$CFLAGS"
        ax_pthread_save_LIBS="$LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"

        # Check for various functions.  We must include pthread.h,
        # since some functions may be macros.  (On the Sequent, we
        # need a special flag -Kthread to make this header compile.)
        # We check for pthread_join because it is in -lpthread on IRIX
        # while pthread_create is in libc.  We check for pthread_attr_init
        # due to DEC craziness with -lpthreads.  We check for
        # pthread_cleanup_push because it is one of the few pthread
        # functions on Solaris that doesn't have a non-functional libc stub.
        # We try pthread_create on general principles.

        AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <pthread.h>
#                       if $ax_pthread_check_cond
#                        error "$ax_pthread_check_macro must be defined"
#                       endif
                        static void routine(void *a) { a = 0; }
                        static void *start_routine(void *a) { return a; }],
                       [pthread_t th; pthread_attr_t attr;
                        pthread_create(&th, 0, start_routine, 0);
                        pthread_join(th, 0);
                        pthread_attr_init(&attr);
                        pthread_cleanup_push(routine, 0);
                        pthread_cleanup_pop(0) /* ; */])],
            [ax_pthread_ok=yes],
            [])

        CFLAGS="$ax_pthread_save_CFLAGS"
        LIBS="$ax_pthread_save_LIBS"

        AC_MSG_RESULT([$ax_pthread_ok])
        AS_IF([test "x$ax_pthread_ok" = "xyes"], [break])

        PTHREAD_LIBS=""
        PTHREAD_CFLAGS=""
done
fi

# Various other checks:
if test "x$ax_pthread_ok" = "xyes"; then
        ax_pthread_save_CFLAGS="$CFLAGS"
        ax_pthread_save_LIBS="$LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"

        # Detect AIX lossage: JOINABLE attribute is called UNDETACHED.
        AC_CACHE_CHECK([for joinable pthread attribute],
            [ax_cv_PTHREAD_JOINABLE_ATTR],
            [ax_cv_PTHREAD_JOINABLE_ATTR=unknown
             for ax_pthread_attr in PTHREAD_CREATE_JOINABLE PTHREAD_CREATE_UNDETACHED; do
                 AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <pthread.h>],
                                                 [int attr = $ax_pthread_attr; return attr /* ; */])],
                                [ax_cv_PTHREAD_JOINABLE_ATTR=$ax_pthread_attr; break],
                                [])
             done
            ])
        AS_IF([test "x$ax_cv_PTHREAD_JOINABLE_ATTR" != "xunknown" && \
               test "x$ax_cv_PTHREAD_JOINABLE_ATTR" != "xPTHREAD_CREATE_JOINABLE" && \
               test "x$ax_pthread_joinable_attr_defined" != "xyes"],
              [AC_DEFINE_UNQUOTED([PTHREAD_CREATE_JOINABLE],
                                  [$ax_cv_PTHREAD_JOINABLE_ATTR],
                                  [Define to necessary symbol if this constant
                                   uses a non-standard name on your system.])
               ax_pthread_joinable_attr_defined=yes
              ])

        AC_CACHE_CHECK([whether more special flags are required for pthreads],
            [ax_cv_PTHREAD_SPECIAL_FLAGS],
            [ax_cv_PTHREAD_SPECIAL_FLAGS=no
             case $host_os in
             solaris*)
             ax_cv_PTHREAD_SPECIAL_FLAGS="-D_POSIX_PTHREAD_SEMANTICS"
             ;;
             esac
            ])
        AS_IF([test "x$ax_cv_PTHREAD_SPECIAL_FLAGS" != "xno" && \
               test "x$ax_pthread_special_flags_added" != "xyes"],
              [PTHREAD_CFLAGS="$ax_cv_PTHREAD_SPECIAL_FLAGS $PTHREAD_CFLAGS"
               ax_pthread_special_flags_added=yes])

        AC_CACHE_CHECK([for PTHREAD_PRIO_INHERIT],
            [ax_cv_PTHREAD_PRIO_INHERIT],
            [AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <pthread.h>]],
                                             [[int i = PTHREAD_PRIO_INHERIT;]])],
                            [ax_cv_PTHREAD_PRIO_INHERIT=yes],
                            [ax_cv_PTHREAD_PRIO_INHERIT=no])
            ])
        AS_IF([test "x$ax_cv_PTHREAD_PRIO_INHERIT" = "xyes" && \
               test "x$ax_pthread_prio_inherit_defined" != "xyes"],
              [AC_DEFINE([HAVE_PTHREAD_PRIO_INHERIT], [1], [Have PTHREAD_PRIO_INHERIT.])
               ax_pthread_prio_inherit_defined=yes
              ])

        CFLAGS="$ax_pthread_save_CFLAGS"
        LIBS="$ax_pthread_save_LIBS"

        # More AIX lossage: compile with *_r variant
        if test "x$GCC" != "xyes"; then
            case $host_os in
                aix*)
                AS_CASE(["x/$CC"],
                    [x*/c89|x*/c89_128|x*/c99|x*/c99_128|x*/cc|x*/cc128|x*/xlc|x*/xlc_v6|x*/xlc128|x*/xlc128_v6],
                    [#handle absolute path differently from PATH based program lookup
                     AS_CASE(["x$CC"],
                         [x/*],
                         [AS_IF([AS_EXECUTABLE_P([${CC}_r])],[PTHREAD_CC="${CC}_r"])],
                         [AC_CHECK_PROGS([PTHREAD_CC],[${CC}_r],[$CC])])])
                ;;
            esac
        fi
fi

test -n "$PTHREAD_CC" || PTHREAD_CC="$CC"

AC_SUBST([PTHREAD_LIBS])
AC_SUBST([PTHREAD_CFLAGS])
AC_SUBST([PTHREAD_CC])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$ax_pthread_ok" = "xyes"; then
        ifelse([$1],,[AC_DEFINE([HAVE_PTHREAD],[1],[Define if you have POSIX threads libraries and header files.])],[$1])
        :
else
        ax_pthread_ok=no
        $2
fi
AC_LANG_POP
])dnl AX_PTHREAD

