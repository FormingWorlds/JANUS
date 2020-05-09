import os
import netCDF4 as net
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import subprocess
import nctools
from subprocess import call
from netCDF4 import Dataset
from atmosphere_column import atmos

matplotlib.use('Agg')  # Or any other X11 back-end

ncfilecff  = net.Dataset('currentlw.cff')
cff   = ncfilecff.variables['cff']

ncfilep  = net.Dataset('profile.p')
p   = ncfilep.variables['p'][:,0,0]

print(np.shape(cff))

for i in range(0,300,10):
    plt.semilogy(cff[i,:,0,0],p/np.max(p))
plt.ylim(1,np.min(p)/np.max(p))
plt.ylabel('p/p$_{s}$')
plt.xlabel('Contribution Function')
plt.savefig('cff-band.pdf')




ncfileuflx  = net.Dataset('currentlw.uflx')
uflx   = ncfileuflx.variables['uflx']

cff_sum = np.sum(cff[:,:,0,0] * uflx[:,-1,0,0][:,None],axis=0)
cff_sum = cff_sum / np.trapz(cff_sum,p/np.max(p),axis=0)


plt.figure()
plt.semilogy(cff_sum,p/np.max(p))
plt.ylim(1,np.min(p)/np.max(p))
plt.ylabel('p/p$_{s}$')
plt.xlabel('Normalised Total Contribution Function')
plt.title('pure_co2.pdf')
plt.savefig('cff-sum.pdf')

