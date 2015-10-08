do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]abs_top_builddir[@],$(abs_top_builddir),g'   \
               -e 's,[@]VERSION[@],$(VERSION),g' \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

EXTRA_DIST = top.in top-build.in
BUILT_SOURCES = top.py top-build

bin_SCRIPTS = top.py top-build

top-build: top-build.in Makefile
	$(do_subst) < $< > $@

top.py: top.in Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
