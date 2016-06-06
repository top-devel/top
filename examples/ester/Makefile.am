do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_SCRIPTS = test.py test-old.py

noinst_DATA = dati eq_ESTER_old_all_lagrange eq_ESTER_all_lagrange d.dat star.h5

EXTRA_DIST = dati.in test.in test-old.in d.in star.in \
			 eq_ESTER_all_lagrange.in eq_ESTER_old_all_lagrange.in
BUILT_SOURCES = dati test.py test-old.py d.dat star.h5 \
				eq_ESTER_all_lagrange eq_ESTER_old_all_lagrange

.in.py: Makefile
	$(do_subst) < $< > $@

.in.dat: Makefile
	cp $< $@

.in.h5: Makefile
	cp $< $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
