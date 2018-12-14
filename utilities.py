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
import sys
import subprocess
from itertools import permutations
from time import sleep
def sciconst():
    scinum = {}
    scinum['e-charge'] = 1.602176565e-19
    scinum['h-plank'] = 6.62607004e-34
    scinum['c-speed'] = 299792458
    return scinum
def ipseed(y,x):
    prime_list = [2]
    v=2
    Jg=True
    while Jg:
        for i in prime_list:
            if v % i == 0:
                 break
            elif i==prime_list[-1]:
                prime_list.append(v)
        v=v+1
        if len(prime_list)==x:
            Jg=False
    prime_list.insert(0,1)
    return prime_list[y:x]
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
def transmission(crystal, SASE):
    scinum =sciconst()
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
def powerz(dct):
    p = dct['power'].astype(np.float)
    p = np.mean(p,axis=0)
    return p
def powers(z0,dct):
    z = dct['z'].astype(np.float)
    n = np.abs(z0-z).argmin()
    p = dct['power'].astype(np.float)
    return p[:,n]
def bunchingz(dct):
    b = dct['bunching'].astype(np.float)
    b = np.mean(b,axis=0)
    return b
def bunchings(z0,dct):
    z = dct['z'].astype(np.float)
    n = np.abs(z0-z).argmin()
    b = dct['bunching'].astype(np.float)
    return b[:,n]
def phases(z0,dct):
    z = dct['z'].astype(np.float)
    n = np.abs(z0-z).argmin()
    b = dct['phi_mid'].astype(np.float)
    return b[:,n]
def Spectrumf(z0, dct):
    scinum = sciconst()
    e_charge = scinum['e-charge']
    h_plank = scinum['h-plank']
    c_speed = scinum['c-speed']
    xlamds = dct['xlamds']
    nslice = dct['nslice']
    zsep = dct['zsep']
    f0 = c_speed/xlamds
    ps = powers(z0,dct)
    phis = phases(z0,dct)
    deltaT = xlamds*zsep/c_speed
    Fs = 1/deltaT                              # Sampling frequency
    Er = 0.01                                  # Energy resolution(ev)
    Fr = e_charge*Er/h_plank                                # Er = h*Fr/e
    Sp = np.round(Fs/Fr)                       # The estimation of sampling point
    if nslice <= Sp:
        P = np.zeros(int(np.round((Sp-nslice)/2)))
        power = np.concatenate((P,ps,P))
        phase = np.concatenate((P,phis,P))
    else:
        power = ps
        phase = phis
    N = len(power)
    Nfft = np.int(2**np.ceil(np.log2(N)))
    Esase = np.sqrt(power)*np.exp(-1j*phase)
    t = np.array(range(len(Esase)))*deltaT
    energy0 = np.trapz(t,abs(Esase)**2)
    Ef = np.fft.fftshift(np.fft.fft(Esase,Nfft))
    f = Fs*(np.array(range(1,Nfft+1))*1.0-(Nfft+1)/2.0)/Nfft+f0
    lamda = c_speed/f
    energy1 = np.trapz(lamda,abs(Ef)**2)
    Nor = abs(energy0/energy1)
    Ef = np.sqrt(Nor)*Ef
    return Ef, f, lamda
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
def current(f2o,step,intype = 'T', flag ='ele'):
    scinum = sciconst()
    c_speed = scinum['c-speed']
    f = open(f2o) 
    alllines = f.readlines()
    f.close()
    m = len(alllines)
    pcharge = float(alllines[2].split()[-1])
    del alllines[0:5]
    m = len(alllines)
    data = ''.join(alllines)
    data = data.split()
    data = np.array(data)
    data = data.reshape(m,5)
    if intype == 'T':
        t = data[:,-1].astype(np.float)
        t = t - np.min(t)
        z = c_speed*t
    if flag == 'ele':
        dt = step/c_speed
        tmp = np.arange(min(t),max(t),dt)
        curI =tmp[0]/m*pcharge/dt
        return z, curI
def moddist(f2o,gamrange,f2w):
    f = open(f2o) 
    alllines = f.readlines()
    f.close()
    headers = alllines[0:5]
    pcharge = headers[2].split()
    psize = headers[3].split()
    del alllines[0:5]
    m = len(alllines)
    data = ''.join(alllines)
    data = data.split()
    data = np.array(data)
    data = data.reshape(m,6)
    gamma = data[:,-1].astype(np.float)
    ind = np.where(gamma > gamrange[0] & gamma < gamrange[1])
    nbeam = alllines[ind]
    psize[-1] = str(len(nbeam))
    psize.append('\n')
    headers[2] = ' '.join(psize)
    pcharge[-1] = '%.6E' %(len(nbeam)/m*float(pcharge[-1]))
    pcharge.append('\n')
    headers[3] = ' '.join(pcharge)
    fo = open(f2w, "w")
    fo.writelines(headers)
    fo.writelines(nbeam)
    fo.close()

def rwwakefield(s,r,s0,tau,rf = 0):
    scinum = sciconst()
    c_speed = scinum['c-speed']
    if any(s<0):
        print('error： s should be >= 0')
        sys.exit()
    Z0 = 120*np.pi
    sig = 2*r^2/Z0/s0^3
    Gamma = c_speed*tau/s0
    if Gamma<=2.5:
        krs0c = np.array([[1.81620118662482, 1.29832152141677],\
                 [0.29540336833708, 0.18173822141741],\
                 [-1.23728482888772, -0.62770698511448],\
                 [1.05410903517018, 0.47383850057072],\
                 [-0.38606826800684, -0.15947258626548],\
                 [0.05234403145974, 0.02034464240245]])
        Qrc = np.array([[1.09524274851589,  2.02903],
               [2.40729067134909,  1.33341005467949],
               [0.06574992723432,  -0.16585375059715],
               [-0.04766884506469,  0.00075548123372]])
        krs0 = 0
        Qr   = 0
        A = [1, np.pi**2/16]
        for j in range(0,6):
            krs0 = krs0 + krs0c[j,rf]*Gamma**j
        for j in range(0,2):
            Qr = Qr + Qrc[j,rf]*Gamma**j
        kr = krs0/s0
        Ew = -A[rf]*Z0*c_speed/np.pi*(np.exp(-kr*s/(2*Qr))*np.cos(kr*s) )/r**2
    else:
        wp = np.sqrt(Z0*c_speed*sig/tau)
        Ew  = -Z0*c_speed/np.pi*(np.exp(-s/(4*c_speed*tau))*np.cos(np.sqrt(2*wp/r/c_speed)*s) )/r**2
    return Ew
def genwake(f2o, shiftDC, f2w):
    scinum = sciconst()
    c_speed = scinum['c-speed']
    currentFile = f2o
    outputFile = f2w
    sig  = 3.5e7  #Al: 'Conductivity (ohm-1*m-1)'
    tau  = 8e-15   #   % Al: relaxation time
    rf   =  1   #              % rf=1: rectangle chamber: rf=0: round chamber    
    r    =2.5 #             % mm, chamber radius 
    Z0 = 120*np.pi
    f = open(currentFile)
    alllines = f.readlines()
    f.close()
    data = ''.join(alllines)
    data = data.split()
    data = np.array(data)
    zs = data[:,0].astype(np.float)
    Ipk = data[:,1].astype(np.float)
    Q = np.trapz(zs/c_speed,Ipk)
    r = r*1E-3
    s0 = (2*r**2/(Z0*sig))**(1/3)
    fI = Ipk/np.trapz(zs,Ipk)
    s = zs - zs[0]
    w = rwwakefield(s,r,s0,tau,rf)
    n = len(s)
    Ew = np.zeros([n,n])
    for j in range(n):
        for i in range(n):
            if i == j:
                break
            else:
                Ew[i,j] = w[j-i]*fI[i]
    dz = np.mean(np.diff(zs))
    Ez = Q*np.sum(Ew,axis=0)*dz
    Ez_mean = np.trapz(zs,fI*Ez)
#    Ez_rms = np.sqrt(np.trapz(zs,fI*(Ez-Ez_mean)**2))
    zs=max(zs)-zs
    zs=np.flipud(zs)
    Ez=np.flipud(Ez)
    Ipk=np.flipud(Ipk)
    if shiftDC:
        Ez = Ez - Ez_mean
    f = open(outputFile)
    headers = ['? VERSION=1.0 \n','? SIZE='+str(len(zs))+' \n','? COLUMNS ZPOS CURPEAK ELOSS \n']
    f.write(headers)
    for i in range(len(zs)):
        line = '%.6f' %(zs[i])+'    '+'%.6f' %(Ipk[i])+'    '+'%.6f' %(Ez[i])+'  \n'
        f.write(line)
    f.close()
def analysis(cores):
    ocmax=312
    bcmax=256
    jg=True
    while jg:

        reto=subprocess.getoutput("bjobs -u all -q beamphysics")
        retb=subprocess.getoutput("bjobs -u all -q beamphysics-mpi")
        og=reto.count('gzhou')
        bg=retb.count('gzhou')
        oa=reto.count('beamphysic')
        ba=retb.count('beamphysic')
        op=reto.count('PEND')
        bp=retb.count('PEND')
#        retg=subprocess.getoutput("bjobs")
#        olan=retg.count('oak001')+retg.count('oak002')
#        blan=retg.count('bullet0001')+retg.count('bullet0002')
#        ocore=4*(retg.count('oak')-olan)
#        bcore=16*(retg.count('bullet')-blan)
        if bp==0:
            print('nobody pending in beamphysics-mpi')
            return 'beamphysics-mpi'
#        if bg < 4:
#            print('gzhou job number is '+str(bg))
#            return 'beamphysics-mpi'
        if op==0:
            print('nobody pending in beamphysics')
            return 'beamphysics'
#        if og==oa and ocore<np.floor((ocmax-2*4)/cores)*cores:
#            print 'beamphysisc-onlyme'
#            return 'beamphysics'
#        if bg==ba and bcore<np.floor((bcmax-2*16)/cores)*cores:
#            print('beamphysics-mpi-onlyme')
#            return 'beamphysics-mpi'
        if og < oa and og*cores < np.floor(ocmax/2):
            print('beamphysics-crowded')
            return 'beamphysics'
        if bg < ba and bg*cores < np.floor(bcmax/2):
            print('beamphysics-mpi-crowded')
            return 'beamphysics-mpi'
        sleep(15)
        print('waiting')
def sub(queue,nodes):
    if queue=='beamphysics-mpi':
        ret=subprocess.getoutput('bsub -a mympi -q beamphysics-mpi -sla bpmpi -o log -n '+str(nodes)+' genesis2_mpi mod.in')
    if queue=='beamphysics':
        ret=subprocess.getoutput('bsub -a mympi -q beamphysics -o log -n '+str(nodes)+' genesis2_mpi mod.in')
    s=ret.split()
    ss=s[1]
    sss=ss[1:-1]
    print(sss+' submitted')
    return sss
def genmodin(ipseed, line, f2o,f2w):
    f=open(f2o)
    alllines=f.readlines()
    f.close()
    alllines[line]= '   IPSEED = '+str(ipseed)+' \n'
    nfile=open(f2w,'w')
    for line in alllines:
        nfile.write(line)
    nfile.close()





    
    
    

    
 
    










































