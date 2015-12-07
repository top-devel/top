do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

dist_noinst_SCRIPTS = test.py test-old.py

dist_noinst_DATA = dati eq_ESTER_old_all_lagrange eq_ESTER_all_lagrange d.dat star.h5

.in.py: Makefile
	$(do_subst) < $< > $@

.in.dat: Makefile
	cp $< $@

.in.h5: Makefile
	cp $< $@

.in: Makefile
	$(do_subst) < $< > $@
