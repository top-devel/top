AM_FFLAGS = -fPIC
AM_FCFLAGS = -fPIC

lib_LTLIBRARIES = libfit.la

libfit_la_SOURCES = curfit.f fpader.f fpback.f fpbspl.f fpchec.f \
					fpcurf.f fpdisc.f fpgivs.f fpknot.f fprati.f fprota.f \
					spalde.f splev.f \
					interpolate.f90 interpolate_derive.f90 \
					interpolate_dderive.f90 smoothen.f90

clean-local:
	rm -f *.mod
