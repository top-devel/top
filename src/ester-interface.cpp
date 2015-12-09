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
        Rs[i] = s.r(-1, i);
    }
}

extern "C"
void get_rr_(int *id, double *rr) {
    cpy_mat(*id, rr, s.r);
}

extern "C"
void get_p_(int *id, double *p) {
    cpy_mat(*id, p, s.p);
}

extern "C"
void get_rho_(int *id, double *rho) {
    cpy_mat(*id, rho, s.rho);
}

extern "C"
void get_g1_(int *id, double *G1) {
    cpy_mat(*id, G1, s.eos.G1);
}

extern "C"
void get_w_(int *id, double *w) {
    cpy_mat(*id, w, s.w);
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
void get_omega_(double *omega) {
    *omega = s.Omega;
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
