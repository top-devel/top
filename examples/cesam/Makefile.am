do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_SCRIPTS = check1D.py test.py check.py

noinst_DATA = dati dati1D \
			  eq_cesam_all_lagrange_avg_1D eq_cesam_all_lagrange_avg \
			  0000-m18a05-ad.osc model-cesam-1D.osc \
			  sismique_ultime_2099couches_tronc.osc order-cesam.inc \
			  order-cesam-1D.inc

EXTRA_DIST = dati.in dati1D.in \
			 eq_cesam_all_lagrange_avg_1D.in eq_cesam_all_lagrange_avg.in \
			 0000-m18a05-ad.osc.in model-cesam-1D.osc.in \
			 sismique_ultime_2099couches_tronc.osc.in order-cesam.inc \
			 test.in check.in check1D.in \
			 order-cesam-1D.inc.in

BUILT_SOURCES = dati dati1D \
				eq_cesam_all_lagrange_avg eq_cesam_all_lagrange_avg_1D \
				0000-m18a05-ad.osc model-cesam-1D.osc \
				sismique_ultime_2099couches_tronc.osc order-cesam.inc \
				check1D.py test.py check.py \
				order-cesam-1D.inc

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
