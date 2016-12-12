#-- 1D lagrangien avg dp=0 --#
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

var dP_P, Er

field pm, g_m, r, rhom, dg_m, rhom_z
scalar Gamma1, Lambda

in

#-------------------------------------------------------------------------------------
equation eqdP_P:
0 = 
    avg(r/Gamma1) * dP_P    +
    avg(r)        * Er'     +
    2             * Er
with (r=1)
    dr(Er, -1) = 0    at r = 0

#-------------------------------------------------------------------------------------
equation eqEr:

0 =
    fp^2                        * Er +
    avg(-pm/rhom)               * dP_P' +
    avg(g_m*(1.0-1.0/Gamma1))   * dP_P +
    avg(-g_m)                   * Er' +
    avg(-dg_m+Lambda*rhom)      * Er
    
with (r=1)
    dr(dP_P, -1) = 0 at r = 1
# this is equivalent to setting delta P/P_0 = 0 on the edge:


