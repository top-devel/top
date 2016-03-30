do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

dist_noinst_SCRIPTS = test.py

dist_noinst_DATA = dati eq_cesam_all_lagrange_avg 0000-m18a05-ad.osc

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@
