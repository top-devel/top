do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_SCRIPTS = test.py check.py

noinst_DATA = dati eq_cesam_all_lagrange_avg 0000-m18a05-ad.osc \
			  sismique_ultime_2099couches_tronc.osc order-cesam.inc

EXTRA_DIST = dati.in eq_cesam_all_lagrange_avg.in 0000-m18a05-ad.osc.in \
			 sismique_ultime_2099couches_tronc.osc.in order-cesam.inc \
			 test.in check.in
BUILT_SOURCES = dati eq_cesam_all_lagrange_avg 0000-m18a05-ad.osc \
				sismique_ultime_2099couches_tronc.osc \
				test.py check.py

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
