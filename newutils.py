# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 10:43:51 2018

@author: GQZhou
"""

#This contains all the functions for wake seed.
#1. read the input of crystall prams 
#2. do transmission
#3. read the SASE fields.
import scipy.io as sio
import numpy as np
#import h5py
import sciconst 
import sys
import subprocess
from itertools import permutations
def readin(f2o):
    f = open(f2o)
    alllines=f.readlines()
    f.close()
    crystal = {}
    for line in alllines:
        if 'thickness' in line:
            crystal['thickness'] = float(line.split('=')[-1])
        if 'bragg' in line:
            crystal['bragg'] = float(line.split('=')[-1])
        if 'asymmetry' in line:
            crystal['asymmetry'] = float(line.split('=')[-1])
        if 'pho_energy' in line:
            crystal['pho_energy'] = float(line.split('=')[-1])
        if 'xr0' in line:
            crystal['xr0'] = float(line.split('=')[-1])
        if 'xi0' in line:
            crystal['xi0'] = float(line.split('=')[-1])
        if 'xrh' in line:
            crystal['xrh'] = float(line.split('=')[-1])
        if 'xih' in line:
            crystal['xih'] = float(line.split('=')[-1])
    crystal['ele_suscept0'] = crystal['xr0'] + 1j*crystal['xi0']
    crystal['ele_susceptH'] = crystal['xrh'] + 1j*crystal['xih']
    crystal['ele_susceptHBar'] = crystal['xrh'] - 1j*crystal['xih']
    return crystal
def readSASE(f2o):
    ftype = f2o.split('.')[-1]
    if ftype == 'mat':
        SASE = sio.loadmat(f2o)
        return SASE
def transmission(crystal, SASE, scinum):
    cry_thickness = crystal['thickness']        
    cry_asymmetry = crystal['asymmetry']     
    ele_suscept0  = crystal['ele_suscept0']         
    ele_susceptH  = crystal['ele_susceptH']  
    ele_susceptHbar  = crystal['ele_susceptHBar']
    h_Plank = scinum['h-plank'] 
    cry_bragg = crystal['bragg']
    pho_energy = crystal['pho_energy']
    e_charge = scinum['e-charge']
    c_speed = scinum['c-speed']
    freq_arry = SASE['freq']
    gamma0 = np.cos(np.deg2rad(cry_bragg+cry_asymmetry-90))          
    gammaH = np.cos(np.deg2rad(cry_bragg-cry_asymmetry+90))          
    asy_factor = gamma0/gammaH                  
    wavelength = h_Plank*c_speed/pho_energy/e_charge 
    ang_freq = 2*np.pi*c_speed/wavelength            
    wave_num = ang_freq/c_speed
    extin_len = np.sqrt(gamma0*np.abs(gammaH))/(wave_num*np.sqrt(ele_susceptH*ele_susceptHbar))
    A = cry_thickness/extin_len
    C = np.exp(1j*ele_suscept0*wave_num*cry_thickness/(2*gamma0))
    G = np.sqrt(np.abs(asy_factor)*(ele_susceptH*ele_susceptHbar))/ele_susceptHbar
    Omega = 2*np.pi*freq_arry-ang_freq
    tmp = -4*Omega*(np.sin(np.deg2rad(cry_bragg))**2/ang_freq)*(1-2*Omega/ang_freq)+ele_suscept0*(1-asy_factor)
    y = wave_num*extin_len/(2*gamma0)*(asy_factor*(tmp))
    Y1 = -y-np.sqrt(y**2+asy_factor/np.abs(asy_factor))
    Y2 = -y+np.sqrt(y**2+asy_factor/np.abs(asy_factor))
    R1 = G*Y1
    R2 = G*Y2
    R00 = np.exp(1j*(ele_suscept0*wave_num*cry_thickness/2/gamma0+A/2*Y1))*(R2-R1)/(R2-R1*np.exp(1j*A/2*(Y1-Y2)))
    R0H = R1*R2*(1-np.exp(1j*A/2*(Y1-Y2)))/(R2-R1*np.exp(1j*A/2*(Y1-Y2)))
    R001 = R00-C
    return R001, R00, R0H
def seed_gen(seed, f2w):
    header1 = '? VERSION=2.0 \n'
    header2 = '? ZPOS          PRADO        ZRAYL       ZWAIST     PHASE \n'
    spower = seed['power']
    sphase = seed['phase']
    zrayl = seed['zrayl']
    zwaist = seed['zwaist']
    zpos = seed['zpos']
    f=open(f2w,'w')
    f.write(header1)
    f.write(header2)
    for i in range(len(spower)):
        s1='%.6E' %(zpos[i])
        s2='%.6E' %(spower[i])
        s3='%.6f' %(zrayl[i])
        s4='%.6f' %(zwaist[i])
        s5='%.6f' %(sphase[i])
        f.write('  '+s1+'   '+s2+'    '+s3+'    '+s4+'    '+s5+'\n')
    f.close()
def miller(a ,phE):
    #in unit of A 3.567
    lambdas = 12.4e3/phE
    nmax = int(np.floor(2*a/lambdas))
    n = range(0,nmax+1)
    m = range(0,nmax+1)
    l = range(0,nmax+1)
    millers = []
    thetas= []
    for i in n:
        for j in m:
            for k in l:
                if lambdas/2/a*np.sqrt(i**2+j**2+k**2)<=1 and lambdas/2/a*np.sqrt(i**2+j**2+k**2) != 0:
                    millers.append([i,j,k])
                    thetas.append(np.arcsin(lambdas/2/a*np.sqrt(i**2+j**2+k**2)))
                #    print(i,j,k)
    cots = 1/np.tan(np.array(thetas))
    #print(cots)
    fmiller = cots.argmin()
    return millers[fmiller], np.rad2deg(thetas[fmiller])
def susceptdb(phE, millerInd, db):
    #millerInd is a list [1,1,1]
    tmp = subprocess.getoutput('ls '+db)
    files = tmp.split('\n')
    for p in permutations(millerInd):
        fname = str(phE)+'eV_'+str(p[0])+'-'+str(p[1])+'-'+str(p[2])
        if fname in files:
            return p
    print('Please insert new data to the DB') 
def ReadGenesis2(fto):
    prm1='    z'
    ftc=open(fto,'r')
    alllines=ftc.readlines()
    ftc.close()
    n=-1
    data=[]
    dct={}
    for line in alllines:
        n += 1
        if 'zsep' in line:
            tmp = line.split('=')[-1]
            tmp1 = tmp.replace('D','E')
            zsep = float(tmp1)            # Read zsep
        if 'xlamds' in line:
            tmp = line.split('=')[-1]
            tmp1 = tmp.replace('D','E')
            xlamds = float(tmp1)        # Read xlamds
        if prm1 in line:
            nn=n
        if 'output: slice' in line:
            break

    dct['xlamds']=xlamds
    dct['zsep']=zsep
    data=alllines[nn+1:n-1]
    m=len(data)
    data=''.join(data)
    data=data.split()
#    data=map(eval,data)
    data=np.array(data)
    data=data.reshape(m,3)
    del alllines[0:n-1]
    title=alllines[6]
    stitle=title.split()
    p=len(stitle)
    nslice=len(alllines)/(m+7)
    dct['nslice']=nslice
    I = []
    for i in range(0,nslice):
        tmp =alllines[m*i+3].split()
        I.append(float(tmp[0]))
        del alllines[m*i:m*i+7]
    I = np.array(I)
    dct['current'] = I
    data1=''.join(alllines)
    del alllines
    data1=data1.split()
  # data1=map(eval,data1)
    data1=np.array(data1)
    if len(data1)!=m*nslice*p:
        print(len(data1),m,nslice,p)
        print('total data numer != nslice*column*row')
        sys.exit()
    data1=data1.reshape(m*nslice,p)
    dct['z']=data[:,0]
    dct['aw']=data[:,1]
    dct['qf']=data[:,2]
    for ii in range(0,len(stitle)):
        dct[stitle[ii]]=data1[:,ii].reshape(nslice,m)
    return dct
def getad(aw,N):#N is the drift length
    mmin = N/(1+aw**2)
    m = np.ceil(mmin)
    ad = np.sqrt(m/N*(1+aw**2)-1)
    return ad
def lclsIIHE_lat(aw_,N_,qf_,slots_,off_):
    class latcell:
        def __init__(self, aw, N, qf):
            self.aw ='           AW         '+'%.6f' %(aw)+'      131.000000        0.000000 \n' \
            +'           AW         0.000000       24.000000        0.000000 \n'
            self.N = N
            self.ad = '           AD         0.000000      131.000000        0.000000 \n' \
            +'           AD         '+'%.6f' %(getad(aw, self.N))+'       24.000000        0.000000 \n'
            
            self.qf = '           QF         0.000000      141.000000        0.000000 \n' \
            +'           QF        '+'%.6f' %(qf)+'        3.000000        0.000000 \n' \
            +'           QF         0.000000       11.000000        0.000000  \n'
    
    class empcell:
        def __init__(self, aw, N, qf):
            self.aw = '           AW         0.000000      '+'%.6f' %(N)+'        0.000000 \n' 
            self.N = N
            self.ad = '           AD         '+'%.6f' %(getad(aw, self.N))+'       '+'%.6f' %(N) \
            +'        0.000000 \n'
            
            self.qf = '           QF         0.000000      141.000000        0.000000 \n' \
            +'           QF        '+'%.6E' %(qf)+'        3.000000        0.000000 \n' \
            +'           QF         0.000000       11.000000        0.000000  \n'
    lat = []
    for i in range(len(aw_)):
        if i in slots_:
            lat.append(empcell(aw_[i],N_[i],qf_[i]))
        else:
            lat.append(latcell(aw_[i],N_[i],qf_[i]))
    del lat[0:off_]
    f  = open('mod.lat','w')
    header1 = '? VERSION = 1.0 \n'
    header2 = '? UNITLENGTH =0.026 \n'
    f.write(header1)
    f.write(header2)
    for i in range(len(lat)):
        f.write(lat[i].aw)
        f.write(lat[i].ad)
        f.write(lat[i].qf)
    f.close()











































