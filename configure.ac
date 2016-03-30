#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
m4_include([version.m4])
AC_INIT([TOP],
        TOP_VERSION)
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.9 no-dist-gzip dist-bzip2 tar-ustar -Wall
                  subdir-objects foreign])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_PREFIX_DEFAULT(${HOME}/local)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])

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

AC_PROG_CXX
AC_PROG_F77($FC)
AC_PROG_LIBTOOL(disable-shared)
AC_PATH_PROG(PERL, perl)

AM_PATH_PYTHON([2.7])
AC_PATH_PROG(F2PY, f2py-$PYTHON_VERSION, f2py)

# Checks for libraries.

AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether Fortran compiler supports procedure interface])
AC_COMPILE_IFELSE(
              [AC_LANG_PROGRAM([], [
               type, abstract :: abstract_type
               contains
               procedure(foo_inter), deferred :: foo
               end type abstract_type
               abstract interface
               subroutine foo_inter(this)
                   import abstract_type
                   class(abstract_type) :: this
               end subroutine foo_inter
               end interface
               ])],
              [inter=yes],
              [inter=no])
AC_MSG_RESULT($inter)
AC_LANG_POP(Fortran)

AS_IF([test "x$inter" == "xno"],
      [AC_ERROR(Your Fortran compiler does not support procedure interface. Try using a newer one)])



AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([ester], [ester_version],
             [found_ester="yes"],
             [found_ester="no"])

AC_LANG_PUSH(C++)
AS_IF([test "x$found_ester" == "xyes"],
      [AC_CHECK_HEADER([ester.h],
                       [LIBS="$LIBS -lester"
                        AC_DEFINE([USE_LIBESTER], [1], [was libester detected])],
                       [found_ester="no"])])
AC_LANG_POP(C++)

AM_CONDITIONAL([USE_ESTER_LIB], [test "x$found_ester" == "xyes"])

AS_IF([test "x$enable_mpi" == "xyes"],
      [AC_SEARCH_LIBS([blacs_setup_], [blacs-openmpi], [],
                      [AC_ERROR(blacs library not found)],
                      [-lblacsF77init-openmpi])
      AC_SEARCH_LIBS([descinit_], [scalapack-openmpi],
                     [LIBS="$LIBS -lblacsF77init-openmpi"],
                     [AC_ERROR(scalapack library not found)])])

AC_CHECK_FUNC(dgemm_     , [have_blas=yes],   [have_blas=no])
AC_CHECK_FUNC(dgetrf_    , [have_lapack=yes], [have_lapack=no])

AS_IF([test "x$have_blas" == "xno" -o "x$have_lapack" == "xno"],
    [AC_SEARCH_LIBS([dgemm_], [mkl_rt],
                    [],
                    [AC_CHECK_LIB([blas], [dgemm_], [],
                                  [AC_ERROR(blas library not found: consider \
                                   setting LDFLAGS)])
                    AC_CHECK_LIB([lapack], [dgetrf_], [],
                                 [AC_ERROR(lapack library not found: consider \
                                  setting LDFLAGS)])],
                                 [-lmkl_intel_thread -lmkl_core -liomp5 -lpthread])])

AC_LANG_PUSH(Fortran)

AC_MSG_CHECKING([whether we can link against the blas library])
AC_LINK_IFELSE(
              [AC_LANG_PROGRAM([], [call dgemm])],
              [link_blas=yes],
              [link_blas=no])
AC_MSG_RESULT($link_blas)

AS_IF([test "x$link_blas" == "xno"],
      [AC_ERROR(Cannot link with blas library. Try setting LIBS...)])

AC_MSG_CHECKING([whether we can link against lapack])
AC_LINK_IFELSE(
              [AC_LANG_PROGRAM([], [call dgetrf])],
              [link_lapack=yes],
              [link_lapack=no])
AC_MSG_RESULT($link_lapack)

AS_IF([test "x$link_lapack" == "xno"],
      [AC_ERROR(Cannot link with lapack library. Try setting BLAS variable)])

AC_LANG_POP(Fortran)

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
                 python/Makefile
                 examples/Makefile
                 examples/ester/Makefile
                 examples/poly_ester/Makefile
                 examples/poly_ester/model/Makefile
                 examples/cesam/Makefile
                 utils/readeq-real-2D-mono
                 utils/readeq-cplx-2D-mono
                 utils/readeq-real-1D-mono
                 utils/readeq-cplx-1D-mono
                 utils/readeq-real-2D-multi
                 utils/readeq-cplx-2D-multi-mpi
                 utils/readeq-real-2D-multi-mpi
                 utils/readeq-cplx-2D-multi
                 src/Makefile])

AC_OUTPUT

cat << EOF

TOP configuration:
------------------
TOP version:    $VERSION
python version: $PYTHON_VERSION
install prefix: $prefix
FC:             $FC
FCFLAGS:        $FCFLAGS
FFLAGS:         $FFLAGS
LDLFAGS:        $LDFLAGS
LIBS:           $LIBS
F2PY:           $F2PY

EOF
