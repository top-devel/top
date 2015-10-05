AM_FCFLAGS = -I. \
			 -I.. \
			 -I@abs_top_builddir@/lib/der \
			 -I@abs_top_builddir@/lib/getpar \
			 -I@abs_top_builddir@/lib/itemlist \
			 -I@abs_top_builddir@/lib/legendre \
			 $(top_flags) \
			 -g

AM_LDFLAGS = -L@prefix@/lib \
			  -Wl,-rpath -Wl,@prefix@/lib

lib_LTLIBRARIES   = libtop.la
libtop_la_LIBADD  = -lgetpar -llegendre -larnoldi -lder -lfit \
					-litemlist -llapack_extra
libtop_la_SOURCES = grid.F90 \
					udirectory.F90 \
					eigensolve.F90 \
					python.F90 \
					hpsort_index.F90 \
					mod_blacs.F90 \
					model.F90 \
					inputs.F90 \
					matrices.F90

bin_PROGRAMS = TOP
TOP_LDADD = -ltop -lgetpar -llegendre -larnoldi -lder -lfit \
			-litemlist -llapack_extra

TOP_SOURCES = TOP.F90 postproc.F90

EXTRA_DIST = order.inc \
			 model_ester.F90 \
			 model_scf.F90 \
			 model_polytrope.F90 \
			 mod_blacs.F90 \
			 hpsort_index.F90 \
			 udirectory.F90 \
			 matrices.F90 \
			 eigensolve.F90 \
			 postproc.F90 \
			 TOP.F90

model.F90:
	cp $(srcdir)/model_$(model).F90 model.F90

model.lo: inputs.lo grid.lo

TOP.o:				postproc.o libtop.la
postproc.lo:		eigensolve.lo udirectory.lo grid.lo
eigensolve.lo:		matrices.lo mod_blacs.lo grid.lo
matrices.lo:		matrices.inc model.lo grid.lo
inputs.lo:			grid.lo grid.lo
udirectory.lo:		inputs.lo model.lo grid.lo
inputs.lo:			mod_blacs.lo grid.lo
python.lo:			inputs.lo eigensolve.lo

clean-local:
	rm -f *.mod *.lo
	rm -f inputs.F90 matrices.inc model.F90 mod_grid.F90

