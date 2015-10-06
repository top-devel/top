do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]PYTHON[@],$(PYTHON),g'


BUILT_SOURCES = top-script.py

.in.py: Makefile createdirtree
	$(do_subst) < $< > $@

.in: Makefile createdirtree
	$(do_subst) < $< > $@

ester_exdir	= @datadir@/top/models/ester
ester_ex_DATA	= models/ester/eq_ESTER_all_lagrange \
				  models/ester/dati \
				  models/ester/d.dat

poly_ester_exdir	= @datadir@/top/models/poly_ester
poly_ester_ex_DATA	= models/poly_ester/eq_poly_ester \
					  models/poly_ester/dati \
					  models/poly_ester/model/domains_boundaries \
					  models/poly_ester/model/enthalpy_et_der \
					  models/poly_ester/model/lambda_eps

scriptdir	= @datadir@/top
script_DATA	= top-script.py

createdirtree:
	mkdir -p $(ester_exdir) $(poly_ester_exdir) \
		models/ester models/poly_ester/model

install-data-hook:
	chmod +x $(DESTDIR)@datadir@/top/top-script.py

EXTRA_DIST = models/ester/eq_ESTER_all_lagrange \
			 models/ester/dati \
			 models/ester/d.dat \
			 models/poly_ester/eq_poly_ester \
			 models/poly_ester/dati \
			 models/poly_ester/model/domains_boundaries \
			 models/poly_ester/model/enthalpy_et_der \
			 models/poly_ester/model/lambda_eps \
			 top-script.in

clean-local:
	rm -rf $(BUILT_SOURCES)
