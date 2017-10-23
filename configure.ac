AC_INIT([asynch], [1.4.2], [samuel-debionne@uiowa.edu])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])

# Check for programms
AX_PROG_CC_MPI([], [], [AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])] )
AC_PROG_CC_STDC #C99
AM_PROG_CC_C_O  #automake < 1.14
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_RANLIB

# Checks for libraries
AX_CHECK_ZLIB
AX_LIB_HDF5([serial])
AX_LIB_POSTGRESQL
AX_BLAS
AX_LAPACK

AC_DEFINE([OS_VERSION], ["Not specified"], [TODO])
AC_DEFINE([CC_VERSION], ["Not specified"], [TODO])
AC_DEFINE([HDF5_VERSION], ["Not specified"], [TODO])
AC_DEFINE([SZ_VERSION], ["Not specified"], [TODO])
AC_DEFINE([ZL_VERSION], ["Not specified"], [TODO])
AC_DEFINE([PQ_VERSION], ["Not specified"], [TODO])

# Checks for packages
PKG_CHECK_MODULES([CHECK], [check >= 0.9.8], [found_libcheck=1], [found_libcheck=0])
PKG_CHECK_MODULES([METIS], [metis >= 5.0.0], [found_metis=1], [found_metis=0])
PKG_CHECK_MODULES([PETSC], [PETSc >= 3.7.0], [found_petsc=1], [found_petsc=0])

# Define automake conditionals
AM_CONDITIONAL([USE_POSTGRESQL], [test "$found_postgresql" = "yes"])
AM_CONDITIONAL([USE_PETSC], [test $found_petsc -eq 1])

# Add preprocessor definition to config
AS_IF([test $found_libcheck -eq 1], [AC_DEFINE([HAVE_LIBCHECK], [1], [Defined if you have libcheck support])])
AS_IF([test $found_metis -eq 1], [AC_DEFINE([HAVE_METIS], [1], [Defined if you have METIS support])])
AS_IF([test $found_petsc -eq 1], [AC_DEFINE([HAVE_PETSC], [1], [Defined if you have PETSc support])])

# Generate config header
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile tools/Makefile])

AC_OUTPUT
