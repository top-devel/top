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
noinst_bin_PROGRAMS = TOP
noinst_bin_SCRIPTS = legpy.so

legpy.so: Makefile legpy.F90 ../lib/legendre/liblegendre.la
	$(F2PY) -c -m legpy						\
		--f90exec=$(FC)						\
		--f77exec=$(F77)					\
		-llegendre 							\
		-L@libdir@ 							\
		-I@abs_top_builddir@/lib/der		\
		-I@abs_top_builddir@/lib/getpar		\
		-I@abs_top_builddir@/lib/itemlist 	\
		-I@abs_top_builddir@/lib/legendre 	\
		--f77flags="@FFLAGS@"  				\
		--f90flags="@FCFLAGS@" 				\
		$(srcdir)/legpy.F90

TOP_LDADD  = -ltop -lgetpar -llegendre -larnoldi -lder -lfit \
				   -litemlist -llapack_extra
TOP_SOURCES = TOP.F90 ester-interface.cpp hpsort_index.F90 grid.F90

model.F90:
	cp $(srcdir)/model_$(model).F90 model.F90

legpy.o:			grid.o

TOP.o:				eigensolve.o mod_blacs.o grid.o
eigensolve.o:		matrices.o mod_blacs.o grid.o
matrices.o:			matrices.inc model.o grid.o config.o
inputs.o:			grid.o
udirectory.o:		inputs.o model.o grid.o
inputs.o:			mod_blacs.o grid.o
model.o:			inputs.o grid.o ester-interface.o abstract_model.o
postproc.o:			eigensolve.o matrices.o grid.o hpsort_index.o

clean-local:
	rm -f *.mod *.lo *.o
	rm -f inputs.F90 matrices.inc model.F90
