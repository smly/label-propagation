#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.64])
AC_PREREQ([2.50])
AC_INIT(lprop, 0.1.3, [kohei-o@is.naist.jp])
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([foreign 1.8.3])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
