AC_INIT([asynch], [1.0], [felipe-quintero@uiowa.edu])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Check for programms
AX_PROG_CC_MPI([], [], [AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])] )
AC_PROG_CC_STDC #C99
AC_PROG_RANLIB

# Checks for libraries
AX_CHECK_ZLIB
AX_LIB_POSTGRESQL
AX_LIB_HDF5([serial])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile forecasters/Makefile tests/Makefile])

# Checks for packages
PKG_CHECK_MODULES([CHECK], [check >= 0.9.8], [HAVE_LIBCHEK=1], [HAVE_LIBCHEK=0])

AC_OUTPUT
