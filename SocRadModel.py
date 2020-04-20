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

def radCompSoc(atm, dirs, recalc):

    # Define path to spectral file
    # spectral_file = dirs["rad_conv"]+"/spectral_files/sp_all_2020/sp_spider"
    spectral_file = dirs["rad_conv"]+"/spectral_files/sp_all_hitran_highres/sp_all_318_hitran"
    # spectral_file = dirs["rad_conv"]+"/spectral_files/sp_all_hitran_lowres/sp_all_318_hitran"

    # print("SOCRATES spectral file:", spectral_file)

    # Solar zenith angle
    zenith_angle    = atm.zenith_angle  # Hamano+15: 54.7, Ranjan+18: 48.2, Katyal+19: 38

    # Surface albedo
    surface_albedo  = atm.albedo_s      # Hamano+15: 0.2, Schaefer+16: 0.75

    # Other parameters
    longitude       = 0
    latitude        = 0
    basis_function  = 1

    # Write values to netcdf
    nctools.ncout_surf('profile.surf', longitude, latitude, basis_function, surface_albedo)
    nctools.ncout2d('profile.tstar', 0, 0, atm.ts, 'tstar', longname="Surface Temperature", units='K')
    nctools.ncout2d('profile.pstar', 0, 0, atm.ps, 'pstar', longname="Surface Pressure", units='PA')
    nctools.ncout2d('profile.szen', 0, 0, zenith_angle, 'szen', longname="Solar zenith angle", units='Degrees')
    nctools.ncout2d('profile.stoa', 0, 0, atm.toa_heating, 'stoa', longname="Solar Irradiance at TOA", units='WM-2')
    # T, P + volatiles
    nctools.ncout3d('profile.t', 0, 0,   atm.p,  atm.tmp, 't', longname="Temperature", units='K')
    nctools.ncout3d('profile.tl', 0, 0,  atm.pl, atm.tmpl, 'tl', longname="Temperature", units='K')
    nctools.ncout3d('profile.p', 0, 0,   atm.p,  atm.p, 'p', longname="Pressure", units='PA')
    nctools.ncout3d('profile.q', 0, 0,   atm.p,  atm.x_gas["H2O"], 'q', longname="q", units='PPMV') 
    nctools.ncout3d('profile.co2', 0, 0, atm.p,  atm.x_gas["CO2"], 'co2', longname="CO2", units='PPMV') 
    nctools.ncout3d('profile.co', 0, 0,  atm.p,  atm.x_gas["CO"], 'co', longname="CO", units='PPMV') 
    nctools.ncout3d('profile.ch4', 0, 0, atm.p,  atm.x_gas["CH4"], 'ch4', longname="CH4", units='PPMV') 
    nctools.ncout3d('profile.h2', 0, 0,  atm.p,  atm.x_gas["H2"], 'h2', longname="H2", units='PPMV') 
    nctools.ncout3d('profile.n2', 0, 0,  atm.p,  atm.x_gas["N2"], 'n2', longname="N2", units='PPMV') 
    nctools.ncout3d('profile.o2', 0, 0,  atm.p,  atm.x_gas["O2"], 'o2', longname="O2", units='PPMV')

    basename = 'profile'
    s = " "

    # Anchor spectral files and run SOCRATES
    seq4 = ("Cl_run_cdf -B", basename,"-s", spectral_file, "-R 1", str(atm.nbands), " -ch ", str(atm.nbands), " -S -g 2 -C 5 -u")
    seq5 = ("fmove", basename,"currentsw")
    seq6 = ("Cl_run_cdf -B", basename,"-s", spectral_file, "-R 1 ", str(atm.nbands), " -ch ", str(atm.nbands), " -I -g 2 -C 5 -u")
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

    # Open netCDF files produced by SOCRATES
    ncfile1  = net.Dataset('currentsw.vflx')
    ncfile2  = net.Dataset('currentsw.sflx')
    ncfile3  = net.Dataset('currentsw.dflx')
    ncfile4  = net.Dataset('currentsw.uflx')
    ncfile5  = net.Dataset('currentsw.nflx')
    ncfile6  = net.Dataset('currentsw.hrts')
    ncfile7  = net.Dataset('currentlw.dflx')
    ncfile8  = net.Dataset('currentlw.nflx')
    ncfile9  = net.Dataset('currentlw.uflx')
    ncfile10 = net.Dataset('currentlw.hrts')

    # Loop through netCDF variables and populate arrays
    vflxsw   = ncfile1.variables['vflx']  # SW downward flux (direct + diffuse)
    uflxsw   = ncfile4.variables['uflx']  # SW upward flux 
    nflxsw   = ncfile5.variables['nflx']  # SW net flux 
    hrtssw   = ncfile6.variables['hrts']  # SW heating rate (K/day)
    dflxlw   = ncfile7.variables['dflx']  # LW downward flux (diffuse)
    uflxlw   = ncfile9.variables['uflx']  # LW upward flux 
    nflxlw   = ncfile8.variables['nflx']  # LW net flux 
    hrtslw   = ncfile10.variables['hrts'] # LW heating rate (K/day)

    ##### Fluxes

    # Upward SW + LW flux summed over all bands (W/m^2)
    atm.SW_flux_up          = np.sum(uflxsw[:,:],axis=0)[:,0,0]
    atm.LW_flux_up          = np.sum(uflxlw[:,:],axis=0)[:,0,0]

    # Downward SW + LW flux summed over all bands (W/m^2)
    atm.SW_flux_down        = np.sum(vflxsw[:,:],axis=0)[:,0,0]
    atm.LW_flux_down        = np.sum(dflxlw[:,:],axis=0)[:,0,0]

    # Net SW + LW flux summed over all bands (W/m^2)
    atm.SW_flux_net         = np.squeeze(np.sum(uflxsw[:,:],axis=0)[:,0,0] - np.sum(vflxsw[:,:],axis=0)[:,0,0])
    atm.LW_flux_net         = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] - np.sum(dflxlw[:,:],axis=0)[:,0,0])

    # Upward SW + LW flux per band, W/m^2/(band)
    atm.LW_spectral_flux_up = uflxlw[:,:,0,0]
    atm.SW_spectral_flux_up = uflxsw[:,:,0,0]

    # Total up- and downward fluxes, (W/m^2)
    atm.flux_up_total       = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] + np.sum(uflxsw[:,:],axis=0)[:,0,0]) 
    atm.flux_down_total     = np.squeeze(np.sum(vflxsw[:,:],axis=0)[:,0,0] + np.sum(dflxlw[:,:],axis=0)[:,0,0])

    # Total net flux (W/m^2)
    atm.net_flux            = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] - np.sum(dflxlw[:,:],axis=0)[:,0,0] + np.sum(uflxsw[:,:],axis=0)[:,0,0] -  np.sum(vflxsw[:,:],axis=0)[:,0,0])

    # Total net flux per band (W/m^2/(band))
    atm.net_spectral_flux   = uflxlw[:,:,0,0] + uflxsw[:,:,0,0] - dflxlw[:,:,0,0] - vflxsw[:,:,0,0]

    ##### Heating rates

    # Heating rates only for no recalc: recalc if tropopause was found
    if recalc == False:

        # Individual heating contributions (K/day)
        atm.SW_heating          = np.sum(hrtssw[:,:],axis=0)[:,0,0]
        atm.LW_heating          = np.sum(hrtslw[:,:],axis=0)[:,0,0]

        # Total heating (K/day)
        atm.net_heating       = np.squeeze(np.sum(hrtssw[:,:],axis=0) + np.sum(hrtslw[:,:],axis=0))

    return atm
