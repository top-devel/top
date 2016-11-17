do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

noinst_SCRIPTS = test.py test-old.py test-non-adia.py check.py check-non-adia.py

noinst_DATA = dati d.dat star.h5 rasalhague.h5 \
			  eq_ESTER_old_all_lagrange eq_ESTER_all_lagrange eq_non_adiabatic \
			  dati-non-adia

EXTRA_DIST = dati.in d.in star.in rasalhague.in \
			 test.in test-old.in test-non-adia.in check.in check-non-adia.py \
			 eq_ESTER_all_lagrange.in eq_ESTER_old_all_lagrange.in \
			 eq_non_adiabatic.in \
			 dati-non-adia
BUILT_SOURCES = dati test.py test-old.py check.py check-non-adia.py d.dat star.h5 \
				eq_ESTER_all_lagrange eq_ESTER_old_all_lagrange \
				eq_non_adiabatic dati-non-adia test-non-adia.py rasalhague.h5

.in.py: Makefile
	$(do_subst) < $< > $@

.in.dat: Makefile
	cp $< $@

.in.h5: Makefile
	cp $< $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
