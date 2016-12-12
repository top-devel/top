#-- 1D euler avg dp=0 --#
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

var Er, dP, Phi, PhiP

field pm, g_m, r, rhom, pm_z, rhom_z, c2
scalar Gamma1, Lambda

in

#-------------------------------------------------------------------------------------
equation eqEr:
0 = 
    avg(r)                  * Er' +
    2.0                     * Er  +
    avg(r/(Gamma1*pm)*pm_z) * Er  +
    avg(r/(Gamma1*pm))      * dP
with (r=1)
    dr(Er, -1) = 0. at r = 0

#-------------------------------------------------------------------------------------
equation eqdP:

dP' =
    fp^2 * avg(rhom)          * Er +
    avg(g_m*(rhom_z-pm_z/c2)) * Er +
    avg(pm_z/(Gamma1*pm))     * dP +
    avg(rhom)                 * PhiP

with (r=1)
    dr(dP, -1) + pm_z*dr(Er,-1) = 0 at r = 1
# this is equivalent to setting delta P/P_0 = 0 on the edge:

#-------------------------------------------------------------------------------------
equation eqPhi:

Phi' = PhiP
with (r=1)
    dr(PhiP, -1) = 0 at r = 0

#-------------------------------------------------------------------------------------
equation eqPhiP:
0 =
    avg(r)                         * PhiP' +
    2.0                            * PhiP  + 
    avg(r*Lambda/c2)               * dP    +
    avg(r*Lambda*(pm_z/c2-rhom_z)) * Er
with (r=1)
    r * dr(PhiP, -1) + dr(Phi, -1) = 0 at r = 1
#-------------------------------------------------------------------------------------
