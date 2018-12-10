#This is used for wake seed generation
import h5py
import numpy as np
import scipy.io as sio
from utils import readin, readSASE, transmission
import matplotlib.pyplot as plt
f2o = 'crystal.pram'
crystal = readin(f2o)
print(crystal)
f2o = 'SASE0.mat'
SASE = readSASE(f2o)
plt.plot(SASE['freq'],SASE['spec'])


