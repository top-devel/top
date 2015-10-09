AM_FCFLAGS = -I. \
			 -I.. \
			 -I@abs_top_builddir@/lib/der \
			 -I@abs_top_builddir@/lib/getpar \
			 -I@abs_top_builddir@/lib/itemlist \
			 -I@abs_top_builddir@/lib/legendre \
			 $(top_flags)

AM_LDFLAGS = -L@libdir@ \
			 -Wl,-rpath -Wl,@libdir@

lib_LTLIBRARIES   = libtop.la
libtop_la_LIBADD  = -lgetpar -llegendre -larnoldi -lder -lfit \
					-litemlist -llapack_extra
libtop_la_SOURCES = grid.F90 \
					eigensolve.F90 \
					python.F90 \
					hpsort_index.F90 \
					mod_blacs.F90 \
					model.F90 \
					inputs.F90 \
					string.F90 \
					matrices.F90

noinst_bindir = @abs_builddir@/debug
noinst_bin_PROGRAMS = top-debug
top_debug_LDADD  = -ltop -lgetpar -llegendre -larnoldi -lder -lfit \
				   -litemlist -llapack_extra
top_debug_SOURCES = TOP.F90

model.F90:
	cp $(srcdir)/model_$(model).F90 model.F90


TOP.o:				libtop.la
postproc.lo:		eigensolve.lo udirectory.lo grid.lo
eigensolve.lo:		matrices.lo mod_blacs.lo grid.lo
matrices.lo:		matrices.inc model.lo grid.lo
inputs.lo:			grid.lo grid.lo
udirectory.lo:		inputs.lo model.lo grid.lo
inputs.lo:			mod_blacs.lo grid.lo string.lo
python.lo:			inputs.lo eigensolve.lo string.lo model.lo
model.lo:			inputs.lo grid.lo

clean-local:
	rm -f *.mod *.lo
	rm -f inputs.F90 matrices.inc model.F90 mod_grid.F90

