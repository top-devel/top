#include "config.h"

#ifdef USE_LIBESTER

#include <ester.h>

star2d s;

extern "C"
int cpp_read_ester_model_(char *model, int *nr, int *nth, int *ndom) {

    s.read(model);
    *nr = s.nr;
    *nth = s.nth;
    *ndom = s.ndomains;

    return 0;
}

void cpy_mat(int id, double *farray, const matrix& field) {
    int skip;
    int idom = id - 1;

    skip = 0;
    for (int i=0; i<idom; i++)
        skip += s.map.gl.npts[i];

    for (int j=0; j<s.nth; j++) {
        for (int i=0; i<s.map.npts[idom]; i++) {
            int idx = j*s.map.gl.npts[idom] + i;
            farray[idx] = field(skip+i, j);
        }
    }
}

extern "C"
void get_rs_(int *id, double *Rs) {
    for (int i=0; i<s.nth; i++) {
        Rs[i] = s.r(-1, i) * s.units.r;
    }
}

extern "C"
void get_rr_(int *id, double *rr) {
    cpy_mat(*id, rr, s.r*s.units.r);
}

extern "C"
void get_p_(int *id, double *p) {
    cpy_mat(*id, p, s.p*s.units.p);
}

extern "C"
void get_rho_(int *id, double *rho) {
    cpy_mat(*id, rho, s.rho*s.units.rho);
}

extern "C"
void get_g1_(int *id, double *G1) {
    cpy_mat(*id, G1, s.eos.G1);
}

extern "C"
void get_g3m_(int *id, double *G3m) {
    cpy_mat(*id, G3m, s.eos.G3_1);
}

extern "C"
void get_cv_(int *id, double *cv) {
    cpy_mat(*id, cv, s.eos.cv);
}

extern "C"
void get_w_(int *id, double *w) {
    cpy_mat(*id, w, s.w*s.units.Omega);
}

extern "C"
void get_kappa_(int *id, double *k) {
    cpy_mat(*id, k, s.opa.k);
}

extern "C"
void get_xi_(int *id, double *xi) {
    cpy_mat(*id, xi, s.opa.xi);
}

extern "C"
void get_dlnxi_lnrho_(int *id, double *dlnxi_lnrho) {
    cpy_mat(*id, dlnxi_lnrho, s.opa.dlnxi_lnrho);
}

extern "C"
void get_dlnxi_lnt_(int *id, double *dlnxi_lnT) {
    cpy_mat(*id, dlnxi_lnT, s.opa.dlnxi_lnT);
}

extern "C"
void get_t_(int *id, double *t) {
    cpy_mat(*id, t, s.T*s.units.T);
}

extern "C"
void get_eps_(int *id, double *eps) {
    cpy_mat(*id, eps, s.nuc.eps);
}

extern "C"
void get_nex_(int *pts) {

    *pts = s.map.nex;
}

extern "C"
void get_npts_(int *pts) {
    for (int i=0; i<s.ndomains; i++) {
        pts[i] = s.map.gl.npts[i];
    }
}

extern "C"
void get_theta_(double *th) {
    for (int i=0; i<s.nth; i++) {
        th[i] = s.th(i);
    }
}

extern "C"
void get_zeta_(double *z) {
    int i;

    for (i=0; i<s.nr; i++) {
        z[i] = s.z(i);
    }
}

extern "C"
void get_mass_(double *m) {
    *m = s.M;
}

extern "C"
void get_radius_(double *r) {
    *r = s.R;
}

extern "C"
void get_lum_(double *lum) {
    *lum = s.luminosity();
}

extern "C"
void get_tc_(double *tc) {
    *tc = s.Tc;
}

extern "C"
void get_omega_(double *omega) {
    *omega = s.Omega * s.units.Omega;
}

extern "C"
void get_x_(double *X) {
    *X = s.X0;
}

extern "C"
void get_z_(double *Z) {
    *Z = s.Z0;
}

extern "C"
void get_xc_(double *Xc) {
    *Xc = s.Xc;
}

extern "C"
void cpp_get_grid_(double *grid, double *th) {
    for(int j=0; j<s.nth; j++) {
        th[j] = s.th(j);
        for(int i=0; i<s.nr; i++) {
            grid[i*s.nth + j] = s.r(i, j);
        }
    }
}

#endif
