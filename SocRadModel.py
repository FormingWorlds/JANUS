'''
SocRadModel.py
Returns heating rates
MDH 25/01/19
'''

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

path_to_socrates = os.getcwd()+"/socrates/socrates_main"

def radCompSoc(atm, toa_heating):

    # Solar zenith angle
    zenith_angle    = 48.2

    # Surface albedo
    surface_albedo  = 0.1

    # Other parameters
    longitude       = 0
    latitude        = 0
    basis_function  = 1

    # Write values to netcdf
    nctools.ncout_surf('profile.surf', longitude, latitude, basis_function, surface_albedo)
    nctools.ncout2d('profile.tstar', 0, 0, atm.ts, 'tstar', longname="Surface Temperature", units='K')
    nctools.ncout2d('profile.pstar', 0, 0, atm.ps, 'pstar', longname="Surface Pressure", units='PA')
    nctools.ncout2d('profile.szen', 0, 0, zenith_angle, 'szen', longname="Solar zenith angle", units='Degrees')
    nctools.ncout2d('profile.stoa', 0, 0, toa_heating, 'stoa', longname="Solar Irradiance at TOA", units='WM-2')
    # T, P, and volatiles
    nctools.ncout3d('profile.t', 0, 0,   atm.p,  atm.tmp, 't', longname="Temperature", units='K')
    nctools.ncout3d('profile.tl', 0, 0,  atm.pl, atm.tmpl, 'tl', longname="Temperature", units='K')
    nctools.ncout3d('profile.p', 0, 0,   atm.p,  atm.p, 'p', longname="Pressure", units='PA')
    nctools.ncout3d('profile.q', 0, 0,   atm.p,  atm.x_gas["H2O"], 'q', longname="q", units='PPMV') 
    nctools.ncout3d('profile.co2', 0, 0, atm.p,  atm.x_gas["CO2"], 'co2', longname="CO2", units='PPMV') 
    nctools.ncout3d('profile.co', 0, 0,  atm.p,  atm.x_gas["CO"], 'co', longname="CO", units='PPMV') 
    nctools.ncout3d('profile.ch4', 0, 0, atm.p,  atm.x_gas["CH4"], 'ch4', longname="ch4", units='PPMV') 
    nctools.ncout3d('profile.h2', 0, 0,  atm.p,  atm.x_gas["H2"], 'h2', longname="H2", units='PPMV') 
    nctools.ncout3d('profile.n2', 0, 0,  atm.p,  atm.x_gas["N2"], 'n2', longname="N2", units='PPMV') 
    nctools.ncout3d('profile.o2', 0, 0,  atm.p,  atm.x_gas["O2"], 'o2', longname="O2", units='PPMV')

    # nctools.ncout3d('profile.h2o',0,0,atm.pl, atm.x_gasl["H2O"],'h2o',longname="h2o",units='PPMV')

    basename = 'profile'
    s = " "

    # Anchor spectral files
    seq4 = ("Cl_run_cdf -B", basename,"-s spectral-files/gen_all_2020/sp_spider -R 1 300 -ch 300 -S -g 2 -C 5 -u")
    seq5 = ("fmove", basename,"currentsw")
    seq6 = ("Cl_run_cdf -B", basename,"-s spectral-files/gen_all_2020/sp_spider -R 1 300 -ch 300 -I -g 2 -C 5 -u")
    seq7 = ("fmove", basename,"currentlw")

    comline1 = s.join(seq4)
    comline2 = s.join(seq5)
    comline3 = s.join(seq6)
    comline4 = s.join(seq7)

    if 1==1:
        os.system(comline1)
        os.system(comline2)
        os.system(comline3)
        os.system(comline4)

    #open netCDF files produced by SOCRATES
    ncfile1 = net.Dataset('currentsw.vflx')
    ncfile2 = net.Dataset('currentsw.sflx')
    ncfile3 = net.Dataset('currentsw.dflx')
    ncfile4 = net.Dataset('currentsw.uflx')
    ncfile5 = net.Dataset('currentsw.nflx')
    ncfile6 = net.Dataset('currentsw.hrts')
    ncfile7 = net.Dataset('currentlw.dflx')
    ncfile8 = net.Dataset('currentlw.nflx')
    ncfile9 = net.Dataset('currentlw.uflx')
    ncfile10 = net.Dataset('currentlw.hrts')

    # Create appropriately sized arrays to hold flux data
    p      = ncfile1.variables['plev'][:]
    levels = len(p)
    vflx   = np.zeros(levels)
    sflx   = np.zeros(levels)
    dflx   = np.zeros(levels)
    uflx   = np.zeros(levels)
    nflx   = np.zeros(levels)
    hrts   = np.zeros(levels-1)
    dflxlw = np.zeros(levels)
    nflxlw = np.zeros(levels)
    uflxlw = np.zeros(levels)
    hrtslw = np.zeros(levels-1)

    # Loop through netCDF variables and populate arrays
    uflxlw = ncfile9.variables['uflx']
    #uflxsw = ncfile4.variables['uflx']
    vflxsw = ncfile1.variables['vflx']
    nflxlw = ncfile8.variables['nflx']
    nflxsw = ncfile5.variables['nflx']
    hrtssw = ncfile6.variables['hrts']
    hrtslw = ncfile10.variables['hrts']

    atm.total_heating = np.flip(np.squeeze(np.sum(hrtssw[:,:],axis=0) + np.sum(hrtslw[:,:],axis=0)))
    # print("Total heating:", atm.total_heating)

    # Sum LW flux over all bands
    atm.LW_flux_up          = np.flip(np.sum(uflxlw[:,:],axis=0)[:,0,0])
    atm.LW_spectral_flux_up = np.flip(uflxlw[:,:,0,0])

    return atm
