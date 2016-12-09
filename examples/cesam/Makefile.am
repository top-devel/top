do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
		   -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_SCRIPTS = check1D.py check-parser-1D.py

noinst_DATA = dati1D \
			  eq_cesam_all_lagrange_avg_1D \
			  0000-m18a05-ad.osc sunmlt-fin-nad.osc \
			  sismique_ultime_2099couches_tronc.osc order-cesam.inc \
			  order-cesam-1D.inc \
			  cesam1D.eq cesam1D_euler_l0.eq

EXTRA_DIST = dati1D.in \
			 eq_cesam_all_lagrange_avg_1D.in \
			 0000-m18a05-ad.osc.in sunmlt-fin-nad.osc.in \
			 sismique_ultime_2099couches_tronc.osc.in order-cesam.inc \
			 test.in check.in check1D.in check-parser-1D.in \
			 order-cesam-1D.inc.in \
			 cesam1D.eq cesam1D_euler_l0.eq

BUILT_SOURCES = dati1D \
				eq_cesam_all_lagrange_avg_1D \
				0000-m18a05-ad.osc sunmlt-fin-nad.osc \
				sismique_ultime_2099couches_tronc.osc order-cesam.inc \
				check1D.py check-parser-1D.py \
				order-cesam-1D.inc \
				cesam1D.eq cesam1D_euler_l0.eq

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
