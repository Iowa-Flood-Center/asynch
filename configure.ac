AC_INIT([asynch], [1.0], [felipe-quintero@uiowa.edu])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Check for programms
AX_PROG_CC_MPI([], [], [AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])] )
AC_PROG_RANLIB

# Checks for libraries
AX_CHECK_ZLIB
AX_LIB_POSTGRESQL

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile forecasters/Makefile])
AC_OUTPUT
