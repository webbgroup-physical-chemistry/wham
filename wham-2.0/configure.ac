#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([wham], 2.0, drew.w.ritchie@gmail.com)
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
CXXFLAGS='-Wall -O3'

# Checks for libraries.
BOOST_REQUIRE([], [
    AC_MSG_ERROR(Boost libararies not counts])
BOOST_PROGRAM_OPTIONS([])

AX_PATH_HDF5_NEW([1.8.0],[yes])

AX_BLAS(,[
    AC_MSG_ERROR(Unable to find blas libraries.)])

AC_CHECK_LIB([m], [main], [], [
    AC_MSG_ERROR(You need the -m math libraries.)])
AC_CHECK_LIB([z], [main], [], [
    AC_MSG_ERROR(You need the -z libgz libraries.)])
#AC_CHECK_LIB([blas], [sgemv], [], [
#    AC_MSG_ERROR(You need the blas libraries.)])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.

AC_CONFIG_FILES([
	src/Makefile
	tools/Makefile
    ])
AC_OUTPUT(Makefile)
