'''
Created 28/01/19

@authors:
Mark Hammond (MH)
Tim Lichtenberg (TL)

SOCRATES radiative-convective model
'''

import numpy as np
import math,phys
import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import matplotlib.pyplot as plt
import matplotlib
import SocRadModel
from atmosphere_column import atmos
import pandas as pd
from scipy import interpolate
import seaborn as sns
import copy

# Constants
L_sun       = 3.828e+26                     # W, IAU definition
AU          = 1.495978707e+11               # m
R_universal = 8.31446261815324            # Universal gas constant, J.K-1.mol-1

# Number of radiation and dry adjustment steps
rad_steps   = 1
conv_steps  = 2

def surf_Planck_nu(atm):
    h   = 6.63e-34
    c   = 3.0e8
    kb  = 1.38e-23
    B   = np.zeros(len(atm.band_centres))
    c1  = 1.191042e-5
    c2  = 1.4387752
    for i in range(len(atm.band_centres)):
        nu      = atm.band_centres[i]
        B[i]    = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))

    B   = B * atm.band_widths/1000.0
    return B

def RadConvEqm(output_dir, time_current, atm, toa_heating, loop_counter, SPIDER_options, standalone):

    atm_dry, atm_moist = radiation_timestepping(atm, toa_heating, rad_steps)

    print("=> Computed OLR (dry, moist):", str(round(atm_dry.LW_flux_up[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2")

    plot_heat_balance(atm_dry, atm_moist)

    # # Write TP and spectral flux profiles for later plotting
    # out_a = np.column_stack( ( atm.temp, atm.p*1e-5 ) ) # K, Pa->bar
    # np.savetxt( output_dir+str(int(time_current))+"_atm_TP_profile.dat", out_a )
    # out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths ) )
    # np.savetxt( output_dir+str(int(time_current))+"_atm_spectral_flux.dat", out_a )

    # Profile used in coupled model
    atm = atm_moist

    return atm.LW_flux_up[0], atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths


# Dry adiabat profile
def dry_adiabat_atm(atm):

    # Calculate cp from molar concentrations
    cp_mix = 0.
    for vol in atm.vol_list.keys():
        cp_mix += atm.vol_list[vol] * ga.cpv(vol)

    # Calculate dry adiabat slope
    atm.Rcp = R_universal / cp_mix

    # Calculate dry adiabat temperature profile for staggered nodes (from ground up)
    for idx, prsl in enumerate(atm.pl):
        atm.tmpl[idx] = atm.ts * ( prsl / atm.ps ) ** ( atm.Rcp )

    # Interpolate temperature from staggered nodes
    atm.tmp = np.interp(atm.p, atm.pl, atm.tmpl)

    return atm
  
# Dry convective adjustment
def DryAdj(atm):

    T   = atm.tmp
    p   = atm.p
    
    # Rcp is global
    # Downward pass
    for i in range(len(T)-1):
        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        
        # Adiabat slope
        pfact = (p1/p2)**atm.Rcp
        
        # If slope is shallower than adiabat (unstable), adjust to adiabat
        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) # Equal layer masses
                              # Not quite compatible with how
                              # heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            atm.tmp[i]   = T1
            atm.tmp[i+1] = T2
    
    # Upward pass
    for i in range(len(T)-2,-1,-1):

        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        pfact = (p1/p2)**atm.Rcp

        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) # Equal layer masses
                              # Not quite compatible with how
                              # heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            atm.tmp[i]   = T1
            atm.tmp[i+1] = T2 

    return atm      

# Moist convective adjustment
def MoistAdj(atm, dT):

    # Apply heating
    tmp_heated = atm.tmp + dT

    # Reset to moist adiabat if convectively unstable
    for idx, tmp_heated in enumerate(tmp_heated):
        atm.tmp[idx] = np.amax([tmp_heated, atm.tmp[idx]])

    return atm 

def plot_heat_balance(atm_dry, atm_moist):

        sns.set_style("ticks")
        sns.despine()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

        ax1.semilogy(atm_dry.tmp,atm_dry.p, color="red", ls="-", label=r'Dry adiabat')
        ax1.semilogy(atm_moist.tmp,atm_moist.p, color="blue", ls="-", label=r'Moist adiabat')

        # ax1.set_xlim([0,np.max(atm.temp)])
        # ax1.set_ylim([np.max(atm.p*1e-5),np.min(atm.p*1e-5)])

        # ax1.set_xticks([0,0.2*np.max(atm.temp),0.4*np.max(atm.temp),0.6*np.max(atm.temp),0.8*np.max(atm.temp),np.max(atm.temp)])

        ax1.legend()
        
        ax2.plot(atm_moist.band_centres,surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray",ls='--',label='Black body ('+str(atm_moist.ts)+" K)")
        ax2.plot(atm_dry.band_centres,atm_dry.LW_spectral_flux_up[:,0]/atm_dry.band_widths, color="red")
        ax2.plot(atm_moist.band_centres,atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths, color="blue")
        ax2.set_xlim([np.min(atm.band_centres),np.max(atm.band_centres)])

        
        # # Wavelength arrays
        # OLR_cm = atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths
        # wavelength  = [ 1e+4*(i**(-1)) for i in atm_moist.band_centres ]    # microns
        # OLR_micron  = [ 1e+4*i for i in OLR_cm ]                            # microns
        # ax2.plot(wavelength, OLR_micron, color="blue")
        # OLR_cm = atm_dry.LW_spectral_flux_up[:,0]/atm_dry.band_widths
        # wavelength  = [ 1e+4*(i**(-1)) for i in atm_dry.band_centres ]      # microns
        # OLR_micron  = [ 1e+4*i for i in OLR_cm ]                            # microns
        # ax2.plot(wavelength, OLR_micron, color="red")


        # # Wave length on top x-axis: http://bit.ly/2Pydvb0
        # def tick_function(wavenumber):
        #     wavelength = [ (1./i)*1e+4 for i in wavenumber ] # to micrometer
        #     return ["%.3f" % z for z in wavelength]
        # ax2b = ax2.twiny()
        # new_tick_locations = [np.min(atm.band_centres), np.min(atm.band_centres)+0.25*(np.max(atm.band_centres)-np.min(atm.band_centres)), np.min(atm.band_centres)+0.5*(np.max(atm.band_centres)-np.min(atm.band_centres)), np.min(atm.band_centres)+0.75*(np.max(atm.band_centres)-np.min(atm.band_centres)), np.max(atm.band_centres)]
        # ax2b.set_xlim(ax2.get_xlim())
        # ax2b.set_xticks(new_tick_locations)
        # ax2b.set_xticklabels(tick_function(new_tick_locations))
        # ax2b.set_xlabel(r"Wavelength ($\mu$m)")

        ax1.invert_yaxis()
        ax1.set_xlabel(r'Temperature $T$ (K)')
        ax1.set_ylabel(r'Pressure $P$ (Pa)')

        # Wavenumber settings
        ax2.set_ylabel(r'OLR (W/m$^2$/cm)')
        ax2.set_xlabel('Wavenumber (1/cm)')
        ax2.legend()
        ax2.set_xlim(left=0, right=4000)
        ax2.set_ylim(bottom=0)

        # # Wavelength settings
        # ax2.set_ylabel(r'OLR (W/m$^2$/$\mu$m)')
        # ax2.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
        # ax2.set_xscale("log")
        # ax2.set_yscale("log") 
        # ax2.set_xlim(right=100)
        # ax2.set_ylim(bottom=1e-14)

        plt.savefig("./output"+'/TP_profile_'+str(round(time_current))+'.pdf', bbox_inches="tight")
        plt.close(fig)
        # print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)))

# Time integration for n steps
def radiation_timestepping(atm, toa_heating, rad_steps):

    # Initialise previous OLR and TOA heating to zero
    PrevOLR_dry         = 0.
    PrevMaxHeat_dry     = 0.
    PrevTemp_dry        = atm.tmp * 0.

    # Build moist adiabat structure
    atm_moist           = ga.general_adiabat(atm)

    # Copy moist pressure arrays for dry adiabat
    atm_dry             = dry_adiabat_atm(copy.deepcopy(atm_moist))

    ### Dry calculation
    # Time stepping
    for i in range(0, rad_steps):

        # Compute radiation, midpoint method time stepping
        atm_dry         = SocRadModel.radCompSoc(atm_dry, toa_heating)
        dT_dry          = atm_dry.total_heating * atm_dry.dt

        # Limit the temperature change per step
        dT_dry          = np.where(dT_dry > 5., 5., dT_dry)
        dT_dry          = np.where(dT_dry < -5., -5., dT_dry)

        # Apply heating
        atm_dry.tmp     += dT_dry

        # # Do the surface balance
        # kturb       = .1
        # atm.tmp[-1] += -atm.dt * kturb * (atm.tmp[-1] - atm.ts)
        
        # Dry convective adjustment
        for iadj in range(conv_steps):
            atm_dry     = DryAdj(atm_dry)

        # Convergence criteria
        dTglobal_dry    = abs(round(np.max(atm_dry.tmp-PrevTemp_dry[:]), 4))
        dTtop_dry       = abs(round(atm_dry.tmp[0]-atm_dry.tmp[1], 4))

        # Inform during runtime
        if i % 2 == 1:
            print("Dry adjustment step", i+1, end=": ")
            print("OLR: " + str(atm_dry.LW_flux_up[0]) + " W/m^2,", "dT_max = " + str(dTglobal_dry) + " K, dT_top = " + str(dTtop_dry) + " K")

        # Reduce timestep if heating is not converging
        if dTglobal_dry < 0.05 or dTtop_dry > 3.0:
            atm_dry.dt  = atm_dry.dt*0.99
            print("Dry adiabat not converging -> dt_new =", atm_dry.dt)

        # Break criteria
        dOLR_dry        = abs(round(atm_dry.LW_flux_up[0]-PrevOLR_dry, 4))
        dbreak_dry      = (0.1*(5.67e-8*atm_dry.ts**4)**0.5)

        # Sensitivity break condition
        if (dOLR_dry < dbreak_dry) and i > 5:
            print("Timestepping break ->", end=" ")
            print("dOLR/step =", dOLR_dry, "W/m^2, dTglobal_dry =", dTglobal_dry)
            break    # break here

        PrevOLR_dry       = atm_dry.LW_flux_up[0]
        PrevMaxHeat_dry   = abs(np.max(atm_dry.total_heating))
        PrevTemp_dry[:]   = atm_dry.tmp[:]

    ### Moist outgoing LW only, no heating

    # Compute radiation, midpoint method time stepping
    atm_moist       = SocRadModel.radCompSoc(atm_moist, toa_heating)

    # dT_moist        = atm_moist.total_heating * atm_moist.dt

    # # Limit the temperature change per step
    # dT_moist        = np.where(dT_moist > 5., 5., dT_moist)
    # dT_moist        = np.where(dT_moist < -5., -5., dT_moist)

    # # Apply heating
    # atm_moist.tmp     += dT_moist

    # # Moist single-step adjustment
    # atm_moist       = MoistAdj(atm_moist, dT_moist)
    
    return atm_dry, atm_moist


def InterpolateStellarLuminosity(star_mass, time_current, time_offset, mean_distance):

    luminosity_df           = pd.read_csv("luminosity_tracks/Lum_m"+str(star_mass)+".txt")
    time_current            = time_current/1e+6         # Myr
    time_offset             = time_offset/1e+6          # Myr
    ages                    = luminosity_df["age"]*1e+3 # Myr
    luminosities            = luminosity_df["lum"]      # L_sol

    # Interpolate luminosity for current time
    interpolate_luminosity  = interpolate.interp1d(ages, luminosities)
    interpolated_luminosity = interpolate_luminosity([time_current+time_offset])*L_sun

    stellar_toa_heating     = interpolated_luminosity / ( 4. * np.pi * (mean_distance*AU)**2. )
    
    return stellar_toa_heating[0], interpolated_luminosity/L_sun

####################################
##### Stand-alone initial conditions
####################################

# Planet age and orbit
time_current  = 1e+7                # yr
time_offset   = 1e+7                # yr
mean_distance = 1.0                 # au

# Surface pressure & temperature
P_surf        = 1e+5              # Pa
T_surf        = 300.               # K

# Volatile molar concentrations: ! must sum to one !
vol_list = { 
              "H2O" : 1.0, 
              "CO2" : .0,
              "H2"  : .0, 
              "N2"  : .0,  
              "CH4" : .0, 
              "O2"  : .0, 
              "CO"  : .0, 
              "He"  : .0,
              "NH3" : .0, 
            }

# Create atmosphere object
atm            = atmos(T_surf, P_surf, vol_list)

# Compute stellar heating
toa_heating, star_luminosity = InterpolateStellarLuminosity(1.0, time_current, time_offset, mean_distance)

# Construct radiative-convective profile + heat flux
LW_flux_up, band_centres, LW_spectral_flux_up_per_band_widths = RadConvEqm("./output", time_current, atm, toa_heating, [], [], standalone=True)
