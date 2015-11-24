AM_FCFLAGS = -I. \
			 -I.. \
			 -fPIC \
			 -I@abs_top_builddir@/lib/der \
			 -I@abs_top_builddir@/lib/getpar \
			 -I@abs_top_builddir@/lib/itemlist \
			 -I@abs_top_builddir@/lib/legendre \
			 $(top_flags)

AM_LDFLAGS = -L@libdir@ \
			 -Wl,-rpath -Wl,@libdir@

# lib_LTLIBRARIES   = libtop.la
# libtop_la_LIBADD  = -lgetpar -llegendre -larnoldi -lder -lfit \
# 					-litemlist -llapack_extra
# libtop_la_SOURCES = grid.F90 \
# 					eigensolve.F90 \
# 					python.F90 \
# 					hpsort_index.F90 \
# 					mod_blacs.F90 \
# 					model.F90 \
# 					inputs.F90 \
# 					string.F90 \
# 					matrices.F90

noinst_bindir = @abs_builddir@/debug
noinst_bin_PROGRAMS = top-debug
top_debug_LDADD  = -ltop -lgetpar -llegendre -larnoldi -lder -lfit \
				   -litemlist -llapack_extra
top_debug_SOURCES = TOP.F90

model.F90:
	cp $(srcdir)/model_$(model).F90 model.F90

TOP.o:				eigensolve.o mod_blacs.o grid.o
# postproc.o:		eigensolve.o udirectory.o grid.o
eigensolve.o:		matrices.o mod_blacs.o grid.o
matrices.o:		matrices.inc model.o grid.o
inputs.o:			grid.o grid.o
udirectory.o:		inputs.o model.o grid.o
inputs.o:			mod_blacs.o grid.o string.o
model.o:			inputs.o grid.o
grid.o:				string.o

clean-local:
	rm -f *.mod *.lo
	rm -f inputs.F90 matrices.inc model.F90 mod_grid.F90

