bin_SCRIPTS = top-build

noinst_SCRIPTS = readeq-real-2D-mono \
				 readeq-cplx-2D-mono \
				 readeq-real-1D-mono \
				 readeq-cplx-1D-mono \
				 readeq-real-2D-multi \
				 readeq-cplx-2D-multi-mpi \
				 readeq-real-2D-multi-mpi \
				 readeq-cplx-2D-multi

CLEANFILES = $(bin_SCRIPTS)
DISTCLEANFILES = $(noinst_SCRIPTS)
EXTRA_DIST = readeq-real-2D-mono \
			 readeq-cplx-2D-mono \
			 readeq-real-1D-mono \
			 readeq-cplx-1D-mono \
			 readeq-real-2D-multi \
			 readeq-cplx-2D-multi-mpi \
			 readeq-real-2D-multi-mpi \
			 readeq-cplx-2D-multi

# if IS_COMPLEX_TYPE
# top_type=1
# else
# top_type=0
# endif
do_subst = sed -e 's,[@]PERL[@],$(PERL),g' \
		   -e 's,[@]IS_COMPLEX_TYPE[@],$(top_type),g'

