from FortranFormat import *
import numpy as np

def read_osc(filename):
    f = open(filename,'r')
    formati = FortranFormat('7I10')
    formatv  = FortranFormat('5D19.12')
    for i in range(6):
        line=f.readline()
    tmp = FortranLine(line, formati)
    (nrmod,iconst,ivar,nbelem)= tmp[0:4]
    const = (nrmod,iconst,ivar,nbelem)
    print '(nrmod,iconst,ivar,nbelem)',(nrmod,iconst,ivar,nbelem)

    glob = np.zeros(iconst)
    var = np.zeros([ivar+nbelem,nrmod])

    nlg=(iconst-1)/5+1 # number of line for "global"

    for il in range(nlg):
        line=f.readline()
        tmp = FortranLine(line, formatv)
        i0 = il*5
        i1 = min(iconst,(il+1)*5)
        glob[i0:i1] = tmp[0:i1-i0]

    print 'glob=',glob

    nlg=(ivar+nbelem-1)/5+1 # number of line for each layer
    for ir in range(nrmod):
        for il in range(nlg):
            line=f.readline()
            tmp = FortranLine(line, formatv)
            i0 = il*5
            i1 = min(ivar+nbelem,(il+1)*5)
            
            var[i0:i1,ir] = tmp[0:i1-i0]

    f.close()
    return const,glob,var

def write_osc(fileout,fileosc,fileatm,const,glob,var):
    formati = FortranFormat('4I10')
    formatv  = FortranFormat('5D19.12')
    (nrmod,iconst,ivar,nbelem)= const

    f = open(fileout,'w')
    f.write('# Temporary model for TOP\n')
    f.write('# patched from '+fileosc+'\n')
    f.write('# and '+fileatm+'\n')
    f.write('# CESAM format\n')
    f.write('0 \n')
    line = FortranLine((nrmod,15,20,0), formati)
    f.write(str(line)+'\n')
    for il in range(3):
        i0 = il*5
        i1 = (il+1)*5
        line = FortranLine(glob[i0:i1], formatv)
        f.write(str(line)+'\n')
    for ir in range(nrmod):
        for il in range(4):
            i0 = il*5
            i1 = (il+1)*5
            line = FortranLine(var[i0:i1,ir], formatv)
            f.write(str(line)+'\n')

    f.close()
    return



fileosc='./top/BUILD/examples/cesam/sismique_ultime_2099couches_tronc.osc'
fileatm='atm.txt'
const,glob,var=read_osc( fileosc )

res=np.loadtxt(fileatm)
h_a=res[:,0]
t_a=res[:,1]
rho_a=res[:,2]
p_a=res[:,3]
pt_a=res[:,4]
cs_a=res[:,5]

write_osc('kk.osc',fileosc,fileatm,const,glob,var)
