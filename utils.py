#This contains all the functions for wake seed.
#1. read the input of crystall prams 
#2. do transmission
#3. read the SASE fields.
import scipy.io as sio
import numpy as np
import h5py
import sciconst 
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
    lambdas = 12.4/phE
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

