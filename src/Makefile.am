AM_FCFLAGS = -I. \
			 -I.. \
			 -I@abs_top_builddir@/lib/der \
			 -I@abs_top_builddir@/lib/getpar \
			 -I@abs_top_builddir@/lib/itemlist \
			 -I@abs_top_builddir@/lib/legendre

AM_LDFLAGS = -lgetpar -llegendre -larnoldi -lder -lfit \
			  -litemlist -llapack_extra \
			  -L@prefix@/lib			\
			  -Wl,-rpath -Wl,@prefix@/lib

TOP: TOP.o \
	mod_grid.o \
	inputs.o \
	model.o \
	udirectory.o \
	matrices.o \
	eigensolve.o \
	postproc.o \
	hpsort_index.o \
	TOP.o \
	mod_blacs.o
	$(FC) $(LIBS) $(LDFLAGS) $(AM_LDFLAGS) $(FCFLAGS) $(AM_FCFLAGS) -o $@ $^

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


.F90.o:
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -c -o $@ $<

model.o: model_$(model).F90
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -c -o $@ $<

TOP.o:			model.o matrices.o inputs.o postproc.o
postproc.o:		eigensolve.o udirectory.o
eigensolve.o:	matrices.o
matrices.o:		matrices.inc model.o
model.o:		inputs.o ../config.h
inputs.o:		mod_grid.o
udirectory.o:	inputs.o
inputs.o:		mod_blacs.o

mod_grid.o: 	mod_grid.F90
mod_blacs.o: 	mod_blacs.F90

clean-local:
	rm -f *.mod
	rm -f inputs.F90 matrices.inc mod_grid.F90

