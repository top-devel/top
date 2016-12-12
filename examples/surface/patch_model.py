from FortranFormat import *
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def read_gong(filename):
    f = open(filename,'r')
    formati = FortranFormat('4I10')
    formatv = FortranFormat('5D16.9')

    for i in range(5):
        line=f.readline()
        print line
    tmp = FortranLine(line, formati)
    (nn, iconst, ivar, ivers)= tmp[0:4]
    const = (nn, iconst, ivar, ivers)
    print '(nn, iconst, ivar, ivers)',const

    glob = np.zeros(iconst)
    var = np.zeros([ivar,nn])
    nlg=(iconst-1)/5+1 # number of line for "global"

    for il in range(nlg):
        line=f.readline()
        tmp = FortranLine(line, formatv)
        i0 = il*5
        i1 = min(iconst,(il+1)*5)
        glob[i0:i1] = tmp[0:i1-i0]

    print 'glob=',glob
    nlg=(ivar-1)/5+1 # number of line for each layer
    for ir in range(nn):
        for il in range(nlg):
            line=f.readline()
            tmp = FortranLine(line, formatv)
            i0 = il*5
            i1 = min(ivar,(il+1)*5)
            
            var[i0:i1,ir] = tmp[0:i1-i0]

    f.close()
    return const,glob,var

def read_atm(fileatm):
    fileavg,fileeos = fileatm
    tab_avg = np.loadtxt(fileavg)
    n_a, nv_a = tab_avg.shape
    tab_eos = np.loadtxt(fileeos)
    n_e, nv_e = tab_eos.shape
    if (n_e != n_a):
        print "ERROR: in files ",fileavg,fileeos,"dimension problem ", n_a, n_e
        return -1
    nv = nv_a + nv_e - 1
    tab = np.zeros((n_a, nv))
    for i in range(n_a):
        tab[i,:nv_a]=tab_avg[i,:]
        tab[i,nv_a:nv_e+nv_a]=tab_eos[i,1:]
    return tab
    
    

def write_osc(fileout,fileosc,fileatm,const,glob,var):
    formati = FortranFormat('4I10')
    formatv = FortranFormat('5D19.12')
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


def patch_model(filegong,fileatm,Teff,filepatch,fileunpatch):

    const,glob,var = read_gong( filegong )
    tab = read_atm(fileatm)

    globo=glob.copy()
    globo[8:10]=0
    varo=var.copy()

    varo[5,:]=0
    varo[15:,:]=0

    (nn, iconst, ivar, ivers)=const

    # write unpatched model within CESAM format
    consto=(nn,iconst,ivar,0)
    write_osc(fileunpatch,filegong,'',consto,globo,varo)

    n_atm = len(tab[:,0])
    n_match = int(n_atm*0.8)
    n_matche = int(n_atm*0.6)

    rho_match = tab[n_match,5]
    r_match = tab[n_match,0]
    r_matche = tab[n_matche,0]
    size_r = r_match-r_matche
    index,=np.where(var[4,:] < rho_match)
    iref = index.max()

    f_r_rho = interp1d(np.log(var[4,iref-3:iref+4]), var[0,iref-3:iref+4], kind='cubic')
    r_shift= f_r_rho(np.log(rho_match))+r_match

    r_atm   = r_shift-tab[:,0]
    r_mod   = var[0,0:iref+2]

    index,=np.where(r_mod > r_mod[iref]+size_r )
    iend=index.max()
    patch = np.linspace(0, 1, iref-iend)

    imod_list=[4,3]
    iatm_list=[5,30]
    var_name=['rho','p gas']

    for iv in range(2):
        imod=imod_list[iv]
        iatm=iatm_list[iv]

        y_atm = tab[:,iatm]
        f_tmp = interp1d(r_atm[::-1], np.log(y_atm[::-1]), kind='cubic')
        varo[imod,0:iref+2]=np.exp(f_tmp(r_mod))
        varo[imod,iend:iref]=np.exp(np.log(var[imod,iend:iref])*patch+np.log(varo[imod,iend:iref])*(1.-patch))

        plt.figure(iv)
        plt.title(var_name[iv])
        plt.semilogy(var[0,:]*1e-5,var[imod,:],'b')
        plt.plot(r_atm*1e-5,tab[:,iatm],'g')
        plt.plot(varo[0,:]*1e-5,varo[imod,:],'r+')
        plt.xlim((r_atm.min()*0.999e-5,r_atm.max()*1.001e-5))
        plt.ylim((tab[:,iatm].min(),tab[:,iatm].max()))
        plt.xlabel('radius [km]')
        plt.legend(['1D','3D','patched'], loc='best')
        plt.show()


    imod_list=[2,9]
    iatm_list=[4,41]
    var_name=['T','Gamma1']

    for iv in range(2):
        imod=imod_list[iv]
        iatm=iatm_list[iv]
        
        y_atm = tab[:,iatm]
        f_tmp = interp1d(r_atm[::-1], y_atm[::-1], kind='cubic')
        varo[imod,0:iref+2]= f_tmp(r_mod)
        varo[imod,iend:iref]=var[imod,iend:iref]*patch+varo[imod,iend:iref]*(1.-patch)

        plt.figure(iv+2)
        plt.title(var_name[iv])
        plt.plot(var[0,:]*1e-5,var[imod,:],'b')
        plt.plot(r_atm*1e-5,tab[:,iatm],'g')
        plt.plot(varo[0,:]*1e-5,varo[imod,:],'r+')
        plt.xlim((r_atm.min()*0.999e-5,r_atm.max()*1.001e-5))
        plt.ylim((tab[:,iatm].min(),tab[:,iatm].max()))
        plt.xlabel('radius [km]')
        plt.legend(['1D','3D','patched'], loc='best')
        plt.show()

    radius_init =var [0,abs(var [2,:]-Teff).argmin()]
    radius_patch=varo[0,abs(varo[2,:]-Teff).argmin()]
    glob[2]=radius_patch
    write_osc(filepatch,filegong,str(fileatm),consto,globo,varo)

    return radius_init, radius_patch


# filegong='./modelS/fgong.l5bi.d.15c'
# fileatm=('./antares/cosc13_tav_table.dat','./antares/cosc13_tav_table_eos.dat')
# Teff = 5777.
# fileunpatch = 'modelS.osc'
# filepatch = 'patch.osc'
# Radu,Radp = patch_model(filegong,fileatm,Teff,filepatch,fileunpatch)
