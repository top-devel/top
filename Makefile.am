if PARSER
PARSER_DIR = parser
else
PARSER_DIR =
endif

SUBDIRS = $(PARSER_DIR) lib utils python examples

EXTRA_DIST = src/Makefile.am src/Makefile.in \
			 src/order-default.inc \
			 src/model_ester.F90 \
			 src/model_scf.F90 \
			 src/model_polytrope.F90 \
			 src/model_poly_ester.F90 \
			 src/model_cesam.F90 \
			 src/mod_blacs.F90 \
			 src/hpsort_index.F90 \
			 src/matrices.F90 \
			 src/eigensolve.F90 \
			 src/python.F90 \
			 src/grid.F90 \
			 src/ester-interface.cpp \
			 src/abstract_model.F90 \
			 src/postproc.F90 \
			 src/legpy.F90 \
			 src/config.F90 \
			 activate-top.sh.in \
			 activate-top.csh.in


docs:
	$(MAKE) html latexpdf -C $(top_srcdir)/doc

test: install
	$(PYTHON) python/top-test

do_subst = sed -e 's,[@]SHELL[@],$(SHELL),g'	\
               -e 's,[@]prefix[@],$(prefix),g'	\
               -e 's,[@]bindir[@],$(bindir),g'	\
               -e 's,[@]libdir[@],$(libdir),g'	\
               -e 's,[@]libexecdir[@],$(libexecdir),g'	\
               -e 's,[@]pyexecdir[@],$(pyexecdir),g'	\
               -e 's,[@]pythondir[@],$(pythondir),g'

noinst_SCRIPTS = activate-top.sh activate-top.csh

activate-top.csh: activate-top.csh.in Makefile
	$(do_subst) < $< > $@

activate-top.sh: activate-top.sh.in Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(noinst_SCRIPTS)

clean-local:
	$(MAKE) clean -C src

distclean-local:
	$(MAKE) distclean -C src
