bin_SCRIPTS = top-ec
noinst_SCRIPTS = readeq

CLEANFILES = $(bin_SCRIPTS)
DISTCLEANFILES = $(noinst_SCRIPTS)
EXTRA_DIST = readeq.in

if IS_COMPLEX_TYPE
top_type=1
else
top_type=0
endif
do_subst = sed -e 's,[@]PERL[@],$(PERL),g' \
		   -e 's,[@]IS_COMPLEX_TYPE[@],$(top_type),g'

readeq: readeq.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@

