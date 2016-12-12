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

var dP, Er

field pm, g_m, r, rhom, pm_z, rhom_z
scalar Gamma1, Lambda

in

#-------------------------------------------------------------------------------------
equation eqdP:
0 = 
    avg(r/(Gamma1*pm)) * dP    +
    avg(r)        * Er'     +
    avg(r*pm_z/(pm*Gamma1)) * Er +
    2             * Er
with (r=1)
    dr(Er, -1) = 0    at r = 0

#-------------------------------------------------------------------------------------
equation eqEr:

0 =
    fp^2                        * Er +
    avg(-1/rhom)                * dP' +
    avg(-g_m/(Gamma1*pm))       * dP +
    avg(g_m*(rhom_z/rhom-pm_z/(Gamma1*pm))+Lambda*rhom)* Er
    
with (r=1)
    dr(dP, -1) + pm_z* dr(Er, -1) = 0 at r = 1
# this is equivalent to setting delta P/P_0 = 0 on the edge:


