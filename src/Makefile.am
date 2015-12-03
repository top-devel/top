AM_FCFLAGS = -I. \
			 -I.. \
			 -fPIC \
			 -I@abs_top_builddir@/lib/der \
			 -I@abs_top_builddir@/lib/getpar \
			 -I@abs_top_builddir@/lib/itemlist \
			 -I@abs_top_builddir@/lib/legendre \
			 $(top_flags)

AM_CXXFLAGS = -fPIC

AM_LDFLAGS = -L@libdir@ \
			 -Wl,-rpath -Wl,@libdir@

noinst_bindir = @abs_builddir@/debug
noinst_bin_PROGRAMS = top-debug

top_debug_LDADD  = -ltop -lgetpar -llegendre -larnoldi -lder -lfit \
				   -litemlist -llapack_extra
top_debug_SOURCES = TOP.F90 ester-interface.cpp hpsort_index.F90

model.F90:
	cp $(srcdir)/model_$(model).F90 model.F90

TOP.o:				eigensolve.o mod_blacs.o grid.o
eigensolve.o:		matrices.o mod_blacs.o grid.o
matrices.o:			matrices.inc model.o grid.o
inputs.o:			grid.o grid.o
udirectory.o:		inputs.o model.o grid.o
inputs.o:			mod_blacs.o grid.o
model.o:			inputs.o grid.o ester-interface.o
model.o:			abstract_model.o
postproc.o:			eigensolve.o matrices.o grid.o hpsort_index.o

clean-local:
	rm -f *.mod *.lo *.o
	rm -f inputs.F90 matrices.inc model.F90
