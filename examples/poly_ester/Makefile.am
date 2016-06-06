SUBDIRS = model

do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_SCRIPTS = test.py check.py

noinst_DATA = dati eq_poly_ester

EXTRA_DIST = dati.in test.in check.in \
			 eq_poly_ester.in
BUILT_SOURCES = dati test.py check.py \
				eq_poly_ester

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
