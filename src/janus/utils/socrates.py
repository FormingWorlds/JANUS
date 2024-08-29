'''
SocRadModel.py
Returns heating rates
MDH 25/01/19
'''

import os, glob, re, shutil
import netCDF4 as net
import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import f90nml

import logging 
log = logging.getLogger("fwl."+__name__)

import janus.utils.nctools as nctools
import janus.utils.RayleighSpectrum as RayleighSpectrum
from janus.utils.atmosphere_column import atmos
import janus.utils.phys as phys

from .. import set_socrates_env  # noqa


def radCompSoc(atm, dirs, recalc, rscatter=False,
               rewrite_cfg=True, rewrite_tmp=True, rewrite_gas=False):
    """Runs SOCRATES to calculate fluxes and heating rates

    Parameters
    ----------
        atm : atmos
            Atmosphere object from atmosphere_column.py
        dirs : dict
            Named directories
        recalc : bool
            Is this function call a 'recalculation' case accounting for a tropopause?
        rscatter : bool
            Include Rayleigh scattering?
        rewrite_cfg : bool
            Re-write configuration values (e.g. TOA heating)
        rewrite_PT : bool
            Re-write temperature and pressure arrays
        rewrite_gas : bool
            Re-write composition files
            
    """

    # Check if bands have been set
    if not atm.bands_set:
        raise Exception("Cannot run radiative transfer because bands have not been loaded into atmos object.")

    # Ask SOCRATES to print info to stdout at runtime
    socrates_print = False

    # Pass namelist information to SOCRATES
    # Only supported from SOCRATES version 2306 onwards (revision 1403)
    socrates_use_namelist = True

    # Define path spectral files
    starspectral_file = dirs["output"]+"star.sf"
    runspectral_file  = dirs["output"]+"runtime.sf"

    # Check that atmosphere is okay
    if np.any(atm.cp <= 0):
        log.error("Negative heat capacity!")
        exit(1)
    vals = []
    for arr in atm.x_gas.values():
        vals.append(arr)
    if np.any(np.array(vals).flatten() < 0):
        log.error("Negative mixing ratio(s)!")
        log.error(str(atm.x_gas))
        exit(1)
    
    # Rayleigh scattering
    if rscatter == True:
        
        # Skip if already exists
        if (not os.path.exists(runspectral_file)) or rewrite_gas:

            shutil.copyfile(starspectral_file,      runspectral_file)
            shutil.copyfile(starspectral_file+"_k", runspectral_file+"_k")

            # Insert Rayleigh scattering into spectral file
            rscatter_allowed= {"CO2", "N2", "H2O"}
            rscatter_snames = []
            rscatter_sratio = []
            for vol in atm.vol_list.keys():
                if vol in rscatter_allowed:
                    vol_lower = str(vol).lower()
                    rscatter_snames.append(vol_lower)
                    rscatter_sratio.append(atm.x_gas[vol][-1])
            RayleighSpectrum.rayleigh_coeff_adder(species_list = rscatter_snames, 
                                                mixing_ratio_list = rscatter_sratio, 
                                                spectral_file_path=runspectral_file,
                                                wavelength_dummy_file_path=dirs["output"]+'wavelength_band_file.txt'
                                                )
        scatter_flag = " -r"
    else:
        scatter_flag = ""

    # If we didn't write a new spectral file above, do it now
    if not os.path.exists(runspectral_file):
        shutil.copyfile(starspectral_file,      runspectral_file)
        shutil.copyfile(starspectral_file+"_k", runspectral_file+"_k")

    # Write values to netcdf: SOCRATES Userguide p. 45
    basename = 'profile'

    check_cfg = lambda fpath: bool(rewrite_cfg or (not os.path.exists(fpath)))
    check_tmp = lambda fpath: bool(rewrite_tmp or (not os.path.exists(fpath)))
    check_gas = lambda fpath: bool(rewrite_gas or (not os.path.exists(fpath)))

    # Write configuration stuff
    fthis = basename+'.surf'
    if check_cfg(fthis): nctools.ncout_surf(fthis, 0, 0, 1, float(atm.albedo_s))

    fthis = basename+'.tstar'
    if check_cfg(fthis) or check_tmp(fthis): nctools.ncout2d(   fthis, 0, 0, atm.ts, 'tstar', longname="Surface Temperature", units='K')

    fthis = basename+'.pstar'
    if check_cfg(fthis): nctools.ncout2d(   fthis, 0, 0, atm.ps, 'pstar', longname="Surface Pressure", units='PA')
    
    fthis = basename+'.szen'
    if check_cfg(fthis): nctools.ncout2d(   fthis, 0, 0, float(atm.zenith_angle), 'szen', longname="Solar zenith angle", units='Degrees')
    
    fthis = basename+'.stoa'
    toah = atm.instellation * (1.0 - atm.albedo_pl) * atm.inst_sf  # Values passed to socrates should not include cos(theta), since it's applied during the RT
    atm.toa_heating = toah * np.cos(atm.zenith_angle * np.pi / 180.0)
    if check_cfg(fthis): nctools.ncout2d(   fthis, 0, 0, toah, 'stoa', longname="Solar Irradiance at TOA", units='WM-2')

    # T, P + volatiles
    fthis = basename+'.t'
    if check_tmp(fthis): nctools.ncout3d(   fthis, 0, 0,   atm.p,  atm.tmp, 't', longname="Temperature", units='K')

    fthis = basename+'.tl'
    if check_tmp(fthis): nctools.ncout3d(   fthis, 0, 0,  atm.pl, atm.tmpl, 'tl', longname="Temperature", units='K')

    fthis = basename+'.p'
    if check_tmp(fthis): nctools.ncout3d(   fthis, 0, 0,   atm.p,  atm.p, 'p', longname="Pressure", units='PA')


    # Gases
    fthis = basename+'.q'
    if check_gas(fthis): nctools.ncout3d(   fthis, 0, 0,   atm.p,  phys.molar_mass['H2O'] / atm.mu * atm.x_gas["H2O"], 'q', longname="q", units='kg/kg') 

    # Clouds
    if atm.do_cloud:
        fthis = basename+'.re'
        if check_tmp(fthis): nctools.ncout3d(   fthis, 0, 0,   atm.p,  atm.re, 're', longname="Effective Radius", units='M') 

        fthis = basename+'.lwm'
        if check_tmp(fthis): nctools.ncout3d(   fthis, 0, 0,   atm.p,  atm.lwm, 'lwm', longname="Liquid Water Mass Fraction", units='GG-1') 

        fthis = basename+'.clfr'
        if check_tmp(fthis): nctools.ncout3d(   fthis, 0, 0,   atm.p,  atm.clfr, 'clfr', longname="Cloud Fraction", units='NONE') 

    allowed_vols = {"CO2", "O3", "N2O", "CO", "CH4", "O2", "NO", "SO2", "NO2", "NH3", "HNO3", "N2", "H2", "He", "OCS"}
    for vol in atm.vol_list.keys():
        if vol in allowed_vols:
            vol_lower = str(vol).lower()
            fthis = basename+'.'+vol_lower
            if check_gas(fthis):
                x_gas_this = phys.molar_mass[vol] / atm.mu * atm.x_gas[vol]
                if np.amax(x_gas_this) > 1.0e-30:
                    nctools.ncout3d(basename+'.'+vol_lower, 0, 0, atm.p,  x_gas_this, vol_lower, longname=vol, units='kg/kg') 

    # Call sequences for run SOCRATES + move data
    if atm.do_cloud:
        seq_sw_ex = ["Cl_run_cdf","-B", basename,"-s", runspectral_file, "-R 1", str(atm.nbands), " -ch ", str(atm.nbands), " -S -g 2 -C ", str(atm.cloud_scheme), " -K ", str(atm.cloud_representation), " -d ", str(atm.droplet_type), " -v ", str(atm.solver), " -u", scatter_flag, " -o"]
        seq_sw_mv = ["fmove", basename,"currentsw"]
        
        seq_lw_ex = ["Cl_run_cdf","-B", basename,"-s", runspectral_file, "-R 1", str(atm.nbands), " -ch ", str(atm.nbands), " -I -g 2 -C ", str(atm.cloud_scheme), " -K ", str(atm.cloud_representation), " -d ", str(atm.droplet_type), " -v ", str(atm.solver), " -u", scatter_flag, " -o"]
        seq_lw_mv = ["fmove", basename,"currentlw"]

    else: # cloud flags require Block 10
        seq_sw_ex = ["Cl_run_cdf","-B", basename,"-s", runspectral_file, "-R 1", str(atm.nbands), " -ch ", str(atm.nbands), " -S -g 2 -C 5 -u", scatter_flag]
        seq_sw_mv = ["fmove", basename,"currentsw"]
        
        seq_lw_ex = ["Cl_run_cdf","-B", basename,"-s", runspectral_file, "-R 1", str(atm.nbands), " -ch ", str(atm.nbands), " -I -g 2 -C 5 -u", scatter_flag]
        seq_lw_mv = ["fmove", basename,"currentlw"]

    # Write namelist file?
    if socrates_use_namelist:

        mmw = 0.0
        for vol in atm.vol_list.keys():
            mmw += ( np.mean(atm.x_gas[vol]) *  phys.molar_mass[vol])   # There's probably a more accurate way to handle this.

        rgas = phys.R_gas/mmw
        cp_avg = np.mean(atm.cp) / mmw

        nml = {
            'socrates_constants': {
                'planet_radius':    atm.planet_radius,      # metres
                'mol_weight_air':   mmw,                    # kg mol-1
                'grav_acc':         atm.grav_s,             # m s-2
                'r_gas_dry':        rgas,                   # J kg-1 K
                'cp_air_dry':       cp_avg                  # J kg-1 K-1
            }
        }

        fthis = basename+".nml"
        if os.path.exists(fthis):
            os.remove(fthis)
        f90nml.write(nml,fthis)
        
        seq_sw_ex.extend(["-N", fthis])
        seq_lw_ex.extend(["-N", fthis])

    # Socrates print to stdout at runtime?
    if socrates_print == True:
        outhandle = None
        seq_sw_ex.extend(["-o"])
        seq_lw_ex.extend(["-o"])
    else:
        outhandle = subprocess.DEVNULL

    # SW calculation with SOCRATES
    subprocess.run(seq_sw_ex,check=True,stderr=outhandle,stdout=outhandle)
    subprocess.run(seq_sw_mv,check=True,stderr=outhandle,stdout=outhandle)

    # LW calculation with SOCRATES
    subprocess.run(seq_lw_ex,check=True,stderr=outhandle,stdout=outhandle)
    subprocess.run(seq_lw_mv,check=True,stderr=outhandle,stdout=outhandle)

    # Open netCDF files produced by SOCRATES
    ncfile1  = net.Dataset('currentsw.vflx')  # direct + diffuse
    ncfile2  = net.Dataset('currentsw.sflx')  # direct
    ncfile3  = net.Dataset('currentsw.dflx')  # diffuse 
    ncfile4  = net.Dataset('currentsw.uflx')  # upward
    ncfile6  = net.Dataset('currentsw.hrts')  # heating rate

    ncfile7  = net.Dataset('currentlw.dflx')  # diffuse
    ncfile9  = net.Dataset('currentlw.uflx')  # upward
    ncfile10 = net.Dataset('currentlw.hrts')  # heating rate

    atm.has_contfunc = os.path.exists('currentlw.cff')
    if atm.has_contfunc:
        ncfile11 = net.Dataset('currentlw.cff')   # contribution function

    # Loop through netCDF variables and populate arrays
    vflxsw   = ncfile1.variables['vflx']  # SW downward flux (direct + diffuse)
    uflxsw   = ncfile4.variables['uflx']  # SW upward flux 
    hrtssw   = ncfile6.variables['hrts']  # SW heating rate (K/day)
    dflxlw   = ncfile7.variables['dflx']  # LW downward flux (diffuse)
    uflxlw   = ncfile9.variables['uflx']  # LW upward flux 
    hrtslw   = ncfile10.variables['hrts'] # LW heating rate (K/day)
    if atm.has_contfunc:
        cff  = ncfile11.variables['cff']  # Contribution function (channel, plev, lat, lon)

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

    # Upward SW + LW flux per band, W/m^2
    atm.LW_spectral_flux_up = uflxlw[:,:,0,0]
    atm.SW_spectral_flux_up = uflxsw[:,:,0,0]

    # Downward SW + LW flux per band, W/m^2
    atm.LW_spectral_flux_down = dflxlw[:,:,0,0]
    atm.SW_spectral_flux_down = vflxsw[:,:,0,0]

    # Total up- and downward fluxes, (W/m^2)
    atm.flux_up_total       = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] + np.sum(uflxsw[:,:],axis=0)[:,0,0]) 
    atm.flux_down_total     = np.squeeze(np.sum(vflxsw[:,:],axis=0)[:,0,0] + np.sum(dflxlw[:,:],axis=0)[:,0,0])

    # Total net flux (W/m^2)
    atm.net_flux            = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] - np.sum(dflxlw[:,:],axis=0)[:,0,0] + np.sum(uflxsw[:,:],axis=0)[:,0,0] -  np.sum(vflxsw[:,:],axis=0)[:,0,0])

    # Total net flux per band (W/m^2)
    atm.net_spectral_flux   = uflxlw[:,:,0,0] + uflxsw[:,:,0,0] - dflxlw[:,:,0,0] - vflxsw[:,:,0,0]
    
    # Contribution function
    if atm.has_contfunc:
        atm.cff   = cff[:,:,0,0]
    # Spectral Upward LW flux
    atm.LW_flux_up_i = uflxlw[:,:,0,0]

        ##### Heating rates

    # Heating rates only for no recalc: recalc if tropopause was found
    if recalc == False:

        # Individual heating contributions (K/day)
        atm.SW_heating          = np.sum(hrtssw[:,:],axis=0)[:,0,0]
        atm.LW_heating          = np.sum(hrtslw[:,:],axis=0)[:,0,0]

        # Total heating (K/day)
        atm.net_heating       = np.squeeze(np.sum(hrtssw[:,:],axis=0) + np.sum(hrtslw[:,:],axis=0))
        atm.net_heating       = np.array(list(atm.net_heating))

    # Close netCDF files
    ncfile1.close()
    ncfile2.close()
    ncfile3.close()
    ncfile4.close()
    ncfile6.close()
    ncfile7.close()
    ncfile9.close()
    ncfile10.close()
    if atm.has_contfunc:
        ncfile11.close()

    return atm

# Sting sorting not based on natsorted package
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

# Clean SOCRATES output files after run
def CleanOutputDir( output_dir ):

    types = ("current*.*", "profile.*") 
    files_to_delete = []
    for files in types:
        files_to_delete.extend(glob.glob(output_dir+"/"+files))

    for file in natural_sort(files_to_delete):
        os.remove(file)

