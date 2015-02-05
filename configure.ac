#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([TOP], [0.8])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.9 no-dist-gzip dist-bzip2 tar-ustar -Wall
                  subdir-objects foreign])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])

# Configure TOP version
AC_ARG_WITH([model],
            AS_HELP_STRING([--with-model],
                           [Select model: ester, scf or polytrope]),
            [AS_IF([test "x$withval" != "xester" -o \
                         "x$withval" != "xscf" -o \
                         "x$withval" != "xpolytrope"],
             [model=$withval],
             [AC_ERROR(--with-model should be ester, scf or polytrope)])],
            [model=ester])
AC_SUBST(model, $model)

AC_ARG_WITH([type],
            AS_HELP_STRING([--with-type],
                           [Select TOP verison: real or complex]),
            [AS_IF([test "x$withval" == "xreal"],
                   [type=real],
                   [AS_IF([test "x$withval" == "xcomplex"],
                          [type=complex],
                          [AC_ERROR(--with-type should be real or complex)])])],
            [type=real])

AS_IF([test "x$type" == "xcomplex"],
      [AC_DEFINE([USE_COMPLEX], [1], [Should the complex version be used])])

AM_CONDITIONAL(IS_COMPLEX_TYPE, [test "x$type" == "xcomplex"])

# MPI support?
AC_ARG_ENABLE([mpi],
              AS_HELP_STRING([--enable-mpi],
                             [Enable MPI support]))
# Checks for programs.
AS_IF([test "x$enable_mpi" == "xyes"],
      [AC_PROG_FC(mpif90 mpiifort)
       AC_DEFINE([USE_MPI], [1], [Define to 1 to enable MPI support])],
      [AC_PROG_FC
       enable_mpi=no])

AC_PROG_F77($FC)
AC_PROG_LIBTOOL
AC_PATH_PROG(PERL, perl)

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AS_IF([test "x$enable_mpi" == "xyes"],
      [AC_SEARCH_LIBS([blacs_setup_], [blacs-openmpi], [],
                      [AC_ERROR(blacs library not found)],
                      [-lblacsF77init-openmpi])
      AC_SEARCH_LIBS([descinit_], [scalapack-openmpi],
                     [LIBS="$LIBS -lblacsF77init-openmpi"],
                     [AC_ERROR(scalapack library not found)])])

AS_IF([test "x${BLAS}" == "x"],
    [AC_SEARCH_LIBS([dgemm_], [mkl_intel_lp64],
                    [LIBS="$LIBS -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"],
                    [AC_CHECK_LIB([blas], [dgemm_], [],
                                  [AC_ERROR(blas library not found: consider \
                                   setting LDFLAGS)])
                    AC_CHECK_LIB([lapack], [dgetrf_], [],
                                 [AC_ERROR(lapack library not found: consider \
                                  setting LDFLAGS)])],
                                 [-lmkl_intel_thread -lmkl_core -liomp5
                                  -lpthread])],
    [LIBS="$LIBS ${BLAS}"])

AC_MSG_CHECKING([blas library])
AC_LINK_IFELSE(
              [AC_LANG_SOURCE([
               void main(void) {dgemm_();}
               ])], [link_blas=yes], [link_blas=no])
AC_MSG_RESULT($link_blas)
AS_IF([test "x$link_blas" == "xno"],
      [AC_ERROR(Cannot link with blas library. Try setting BLAS variable)])
AC_MSG_CHECKING([blas library])

AC_LINK_IFELSE(
              [AC_LANG_SOURCE([
               void main(void) {dgetrf_();}
               ])], [link_lapack=yes], [link_lapack=no])
AC_MSG_RESULT($link_lapack)
AS_IF([test "x$link_lapack" == "xno"],
      [AC_ERROR(Cannot link with lapack library. Try setting BLAS variable)])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 lib/Makefile
                 lib/getpar/Makefile
                 lib/legendre/Makefile
                 lib/der/Makefile
                 lib/arnoldi/Makefile
                 lib/fit/Makefile
                 lib/itemlist/Makefile
                 lib/lapack_extra/Makefile
                 utils/Makefile
                 utils/top-ec
                 src/Makefile])
AC_OUTPUT

cat << EOF

--------------------------------------------------------------------------------
TOP build configuration:

FC:         $FC
FCFLAGS:    $FCFLAGS
LIBS:       $LIBS
prefix:     $prefix

TOP version:
model:      $model
type:       $type
dimension:  2D
parallel:   $enable_mpi
--------------------------------------------------------------------------------

EOF
