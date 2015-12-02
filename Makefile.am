SUBDIRS = lib utils python examples

EXTRA_DIST = src/Makefile.am src/Makefile.in \
			 src/order.inc \
			 src/model_ester.F90 \
			 src/model_scf.F90 \
			 src/model_polytrope.F90 \
			 src/model_poly_ester.F90 \
			 src/mod_blacs.F90 \
			 src/hpsort_index.F90 \
			 src/matrices.F90 \
			 src/eigensolve.F90 \
			 src/udirectory.F90 \
			 src/python.F90 \
			 src/grid.F90 \
			 src/string.F90


DISTCLEANFILES = src/Makefile

docs:
	$(MAKE) html latexpdf -C $(top_srcdir)/doc

do_subst = sed -e 's,[@]SHELL[@],$(SHELL),g'	\
               -e 's,[@]prefix[@],$(prefix),g'	\
               -e 's,[@]bindir[@],$(bindir),g'	\
               -e 's,[@]libdir[@],$(libdir),g'	\
               -e 's,[@]libexecdir[@],$(libexecdir),g'	\
               -e 's,[@]pyexecdir[@],$(pyexecdir),g'	\
               -e 's,[@]pythondir[@],$(pythondir),g'

noinst_SCRIPTS = activate-top.sh

activate-top.sh: sourceme.in Makefile
	$(do_subst) < $< > $@
