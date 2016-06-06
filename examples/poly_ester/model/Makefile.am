do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_DATA = domains_boundaries enthalpy_et_der lambda_eps
EXTRA_DIST = domains_boundaries.in enthalpy_et_der.in lambda_eps.in
BUILT_SOURCES = domains_boundaries enthalpy_et_der lambda_eps

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
