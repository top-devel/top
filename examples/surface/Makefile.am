do_subst = sed -e 's,[@]builddir[@],$(abs_builddir),g' \
		   -e 's,[@]PYTHON[@],$(PYTHON),g' \
		   -e 's,[@]top_builddir[@],$(abs_top_builddir),g' \
		   -e 's,[@]srcdir[@],$(abs_srcdir),g'

noinst_SCRIPTS = script.sh run_model.py

noinst_DATA = cesam1D_eul_l0.eq cesam1D_eul_iso_l0.eq \
			  order-cesam-1D_l0.inc order-cesam-1D_l0.inc

EXTRA_DIST = script.in run_model.in \
			 cesam1D_eul_l0.eq.in cesam1D_eul_iso_l0.eq.in \
			 order-cesam-1D_l0.inc.in order-cesam-1D_l0.inc.in


BUILT_SOURCES = script.sh run_model.py \
			 cesam1D_eul_l0.eq cesam1D_eul_iso_l0.eq \
			 order-cesam-1D_l0.inc order-cesam-1D_l0.inc.in


.in.sh: Makefile
	$(do_subst) < $< > $@
	chmod +x $@

.in.py: Makefile
	$(do_subst) < $< > $@

.in: Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
