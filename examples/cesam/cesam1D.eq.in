input double mass
input double Rota
input double age
input string modelfile
input string gridfile
input int pert_model
input string grid_type
input double C0
input double C1
input double C2
input double C3

var dP_P, (Er, Et), PhiP, Phi

field pm, g_m, r, rhom, dg_m, rhom_z
scalar Gamma1, Lambda

in

#-------------------------------------------------------------------------------------
equation eqdP_P:
lh*(lh+1) * Et = 
    avg(r/Gamma1) * dP_P    +
    avg(r)        * Er'     +
    2             * Er
with (r=1)
    dr(Er, -1) = lh * dr(Et, -1)
    at r = 0

#-------------------------------------------------------------------------------------
equation eqEr:

PhiP =
    fp^2                        * Er +
    avg(-pm/rhom)               * dP_P' +
    avg(g_m*(1.0-1.0/Gamma1))   * dP_P +
    avg(-g_m)                   * Er' +
    avg(-dg_m)                  * Er
with (r=1)
    dr(dP_P, -1) = 0 at r = 1
# this is equivalent to setting delta P/P_0 = 0 on the edge:
# termbc s nr w0 1d0                      dP_P^-1(nr)

#-------------------------------------------------------------------------------------
equation eqEt:

0 =
    fp^2 * r    * dr(Et, -1) +
    # 1    * dr(Et, -1) +
# instruction if (lh.eq.0) $prev = 1d0
    -pm/rhom    * dr(dP_P, -1) -
    1.0         * dr(Phi, -1) -
    g_m         * dr(Er, -1)

#-------------------------------------------------------------------------------------
equation eqPhi:

PhiP = Phi'
with (r=1)
    lh * dr(Phi, -1) + (-r + max(0,1-lh)) * dr(PhiP, -1) = 0 at r = 0
#-------------------------------------------------------------------------------------
equation eqPhiP:

0 = 
    avg(r^2)                        * PhiP' +
    avg(2.0*r)                      * PhiP -
    (lh*(lh+1))                     * Phi +
    avg(-r^2*Lambda*rhom/Gamma1)    * dP_P +
    avg(r^2*Lambda*rhom_z)          * Er
with (r=1)
    r * dr(PhiP, -1) + dble(lh+1) * dr(Phi, -1) = 0 at r = 1
# #-------------------------------------------------------------------------------------
