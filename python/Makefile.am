do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

bin_SCRIPTS = top.py

BUILT_SOURCES = top.py

.in.py: Makefile
	$(do_subst) < $< > $@
