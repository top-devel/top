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

void cpy_mat(int id, double *farray, matrix& field) {
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
void get_rho_(int *id, double *rho) {
    cpy_mat(*id, rho, s.rho);
}

extern "C"
void get_nex_(int *pts) {

    *pts = s.map.nex;
}

extern "C"
void get_npts_(int *pts) {
    int i;

    for (i=0; i<s.ndomains; i++) {
        pts[i] = s.map.gl.npts[i];
    }
}

extern "C"
void get_theta_(double *th) {
    int i;

    for (i=0; i<s.nth; i++) {
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
