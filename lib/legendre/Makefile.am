lib_LTLIBRARIES = liblegendre.la

liblegendre_la_SOURCES = gauleg.f90 spherics.f90 fast_legendre.f90 \
						 legendre.f90 fast_pylm.f90 integrales.f90

legendre.lo:	fast_legendre.lo
integrales.lo:	fast_pylm.lo
spherics.lo:	fast_legendre.lo fast_pylm.lo

clean-local:
	rm -f *.mod
