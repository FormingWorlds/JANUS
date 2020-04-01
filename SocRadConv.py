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
conv_steps  = 0

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
    B   = (1.-atm.albedo_s) * np.pi * B * atm.band_widths/1000.0
    return B

# def SpectralExitance_nu(atm):
#     h   = 6.62607015e-34
#     c   = 2.99792458e+8
#     kb  = 1.380649e-23
#     B   = np.zeros(len(atm.band_centres))
#     for i, nu in enumerate(atm.band_centres[i]):
#         B[i] = ( 2.*h*(c**2.)*(nu**3) / ( np.exp( h*c*nu / (kb*atm.ts) ) - 1 ) )
#     B   *= 
#     B   *= atm.band_widths
#     return B

def RadConvEqm(output_dir, star_age, atm, toa_heating, loop_counter, SPIDER_options, standalone, cp_dry):

    atm_dry, atm_moist = radiation_timestepping(atm, toa_heating, rad_steps, cp_dry)

    # Inform user + plot
    if standalone == True:
        plot_flux_balance(atm_dry, atm_moist, cp_dry, star_age)
        print("Computed OLR => moist:", str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2", end=" ")
        if cp_dry == True:
            print("| dry:", str(round(atm_dry.LW_flux_up[0], 3)) + " W/m^2", end=" ")
        print()

    # # Write TP and spectral flux profiles for later plotting
    # out_a = np.column_stack( ( atm.temp, atm.p*1e-5 ) ) # K, Pa->bar
    # np.savetxt( output_dir+str(int(time_current))+"_atm_TP_profile.dat", out_a )
    # out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths ) )
    # np.savetxt( output_dir+str(int(time_current))+"_atm_spectral_flux.dat", out_a )

    # Profile used in coupled model
    atm = atm_moist

    # return atm.LW_flux_up[0], atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths
    return atm

# Dry adiabat profile
def dry_adiabat_atm(atm):

    # Calculate cp from molar concentrations
    cp_mix = 0.
    for vol in atm.vol_list.keys():
        cp_mix += atm.vol_list[vol] * ga.cpv(vol, atm.ts)

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

def plot_flux_balance(atm_dry, atm_moist, cp_dry, star_age):

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13,10))
        sns.set_style("ticks")
        sns.despine()

        # Line settings
        col_idx  = 3
        col_vol1 = "H2O"
        col_vol2 = "N2"
        col_vol3 = "H2"
        col_vol4 = "O2"

        # Temperature vs. pressure
        ax1.semilogy(atm_moist.tmp,atm_moist.p, color=ga.vol_colors[col_vol1][col_idx+1], ls="-", label=r'Moist adiabat')
        if cp_dry == True: ax1.semilogy(atm_dry.tmp,atm_dry.p, color=ga.vol_colors[col_vol3][col_idx+1], ls="-", label=r'Dry adiabat')
        ax1.legend()
        ax1.invert_yaxis()
        ax1.set_xlabel(r'Temperature $T$ (K)')
        ax1.set_ylabel(r'Pressure $P$ (Pa)')
        ax1.set_ylim(bottom=atm_moist.ps*1.01)

        # Fluxes vs. pressure

        # Zero line
        ax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
        ax2.axvline(-1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)
        ax2.axvline(1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)

        # LW down
        if cp_dry == True: ax2.semilogy(atm_dry.LW_flux_down*(-1),atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx+1], ls=(0, (3, 1, 1, 1)))
        ax2.semilogy(atm_moist.LW_flux_down*(-1),atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx+1], ls=(0, (3, 1, 1, 1)), label=r'$F_\mathrm{LW}^{\downarrow}$')
        
        # SW down
        if cp_dry == True: ax2.semilogy(atm_dry.SW_flux_down*(-1),atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx-1], ls="--")
        ax2.semilogy(atm_moist.SW_flux_down*(-1),atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx-1], ls="--", label=r'$F_\mathrm{SW}^{\downarrow}$')
        
        # Net flux
        if cp_dry == True: ax2.semilogy(atm_dry.net_flux,atm_dry.pl, color=ga.vol_colors[col_vol3][6], ls="-", lw=2)
        ax2.semilogy(atm_moist.net_flux,atm_moist.pl, color=ga.vol_colors[col_vol1][6], ls="-", lw=2, label=r'$F_\mathrm{net}$')
    
        # LW_up - SW_down
        if cp_dry == True: ax2.semilogy(atm_dry.LW_flux_up-atm_dry.SW_flux_down,atm_dry.pl, color=ga.vol_colors[col_vol4][col_idx], ls="-")
        net_heating_profile = atm_moist.LW_flux_up - atm_moist.SW_flux_down
        ax2.semilogy(net_heating_profile,atm_moist.pl, color=ga.vol_colors[col_vol2][col_idx], ls="-", label=r'$F_\mathrm{LW}^{\uparrow}$-$F_\mathrm{SW}^{\downarrow}$')

        # Plot tropopause
        trpp_idx = int(atm_moist.trpp[0])
        if trpp_idx > 0:
            ax2.axhline(atm_moist.pl[trpp_idx], color=ga.vol_colors[col_vol2][col_idx], lw=0.5, ls="--")

        # SW up
        if cp_dry == True: ax2.semilogy(atm_dry.SW_flux_up,atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx-1], ls=":")
        ax2.semilogy(atm_moist.SW_flux_up,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx-1], ls=":", label=r'$F_\mathrm{SW}^{\uparrow}$')
        
        # LW up
        if cp_dry == True: ax2.semilogy(atm_dry.LW_flux_up,atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx+1], ls=(0, (5, 1)))
        ax2.semilogy(atm_moist.LW_flux_up,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx+1], ls=(0, (5, 1)), label=r'$F_\mathrm{LW}^{\uparrow}$')

        # Heating rate
        # ax2.semilogy(atm_moist.total_heating,atm_moist.p, color=ga.vol_colors[col_vol2][1], ls=":", label=r'Total heating')
        
        ax2.legend(ncol=6, fontsize=8, loc=3)
        ax2.invert_yaxis()
        ax2.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
        ax2.set_xlabel(r'Outgoing flux $F^{\uparrow}$ (W m$^{-2}$)')
        ax2.set_ylabel(r'Pressure $P$ (Pa)')
        ax2.set_ylim(bottom=atm_moist.ps*1.01)

        # Wavenumber vs. OLR
        ax3.plot(atm_moist.band_centres,surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray",ls='--',label=str(round(atm_moist.ts))+' K blackbody')
        if cp_dry == True: ax3.plot(atm_dry.band_centres,atm_dry.LW_spectral_flux_up[:,0]/atm_dry.band_widths, color=ga.vol_colors[col_vol2][col_idx+1])
        ax3.plot(atm_moist.band_centres,atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol1][col_idx+1])
        ax3.set_xlim([np.min(atm.band_centres),np.max(atm.band_centres)])
        ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)')
        ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
        ax3.legend()
        ax3.set_xlim(left=0, right=5000)
        ax3.set_ylim(bottom=0)
        
        # Wavelength versus OLR log plot
        OLR_cm_moist = atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths
        wavelength_moist  = [ 1e+4/i for i in atm_moist.band_centres ]          # microns
        OLR_micron_moist  = [ 1e+4*i for i in OLR_cm_moist ]                    # microns
        if cp_dry == True: 
            OLR_cm_dry = atm_dry.LW_spectral_flux_up[:,0]/atm_dry.band_widths
            wavelength_dry  = [ 1e+4/i for i in atm_dry.band_centres ]              # microns
            OLR_micron_dry  = [ 1e+4*i for i in OLR_cm_dry ]                        # microns
            ax4.plot(wavelength_dry, OLR_micron_dry, color=ga.vol_colors[col_vol3][col_idx+1])
        
        ax4.plot(wavelength_moist, OLR_micron_moist, color=ga.vol_colors[col_vol1][col_idx+1])
        ax4.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)')
        ax4.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
        ax4.set_xscale("log")
        ax4.set_yscale("log") 
        ax4.set_xlim(left=0.1, right=100)
        ax4.set_ylim(bottom=1e-20, top=1e5)
        # ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])
        ax4.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100])
        ax4.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])

        plt.savefig("./output"+'/TP_profile_'+str(round(star_age))+'.pdf', bbox_inches="tight")
        plt.close(fig)

# Time integration for n steps
def radiation_timestepping(atm, toa_heating, rad_steps, cp_dry):

    # Initialise previous OLR and TOA heating to zero
    PrevOLR_dry         = 0.
    PrevMaxHeat_dry     = 0.
    PrevTemp_dry        = atm.tmp * 0.

    # Build moist adiabat structure
    atm_moist           = ga.general_adiabat(atm)



    # Copy moist pressure arrays for dry adiabat
    atm_dry             = dry_adiabat_atm(copy.deepcopy(atm_moist))

    if cp_dry == True:

        ### Dry calculation
        # Time stepping
        for i in range(0, rad_steps):

            # toa_heating = 0

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
                print("Dry adiabat not converging -> dt_new =", atm_dry.dt, "days")

            # Break criteria
            dOLR_dry        = abs(round(atm_dry.LW_flux_up[0]-PrevOLR_dry, 6))
            dbreak_dry      = (0.01*(5.67e-8*atm_dry.ts**4)**0.5)

            # Sensitivity break condition
            if (dOLR_dry < dbreak_dry) and i > 5:
                print("Timestepping break ->", end=" ")
                print("dOLR/step =", dOLR_dry, "W/m^2, dTglobal_dry =", dTglobal_dry)
                break    # break here

            PrevOLR_dry       = atm_dry.LW_flux_up[0]
            PrevMaxHeat_dry   = abs(np.max(atm_dry.total_heating))
            PrevTemp_dry[:]   = atm_dry.tmp[:]

    ### Moist outgoing LW only, no heating

    # Compute radiation, no moist timestepping
    atm_moist       = SocRadModel.radCompSoc(atm_moist, toa_heating)

    print("OLR w/o stratosphere:", str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    # Redo calculation w/ stratosphere
    net_LW_SW = atm_moist.LW_flux_up - atm_moist.SW_flux_down

    # Find tropopause index
    signchange = ((np.roll(np.sign(net_LW_SW), 1) - np.sign(net_LW_SW)) != 0).astype(int)
    signchange[0] = 0
    
    if len(np.nonzero(signchange)) > 1:
        print("Error: found more than one tropopause solution!")
    elif (np.max(signchange) == 0):
        print("No tropopause")
    
    # Adjust all values above tropopause level    
    else:

        # Tropopause index, pressure and temperature
        trpp_idx          = np.nonzero(signchange)[0][0]

        # Inform user
        print("Tropopause @ (index, P/Pa, T/K):", trpp_idx, round(atm_moist.pl[trpp_idx],2), round(atm_moist.tmpl[trpp_idx],2))

        atm_moist.trpp[0] = trpp_idx                  # index
        atm_moist.trpp[1] = atm_moist.pl[trpp_idx]    # pressure 
        atm_moist.trpp[2] = atm_moist.tmpl[trpp_idx]  # temperature
        
        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        # Rerun SOCRATES on new atmosphere structure
        atm_moist       = SocRadModel.radCompSoc(atm_moist, toa_heating)

        print("OLR w/ stratosphere:", str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_dry, atm_moist

def set_stratosphere(atm):

    trpp_idx = int(atm.trpp[0])
    trpp_prs = atm.trpp[1]
    trpp_tmp = atm.trpp[2]

    # Standard nodes
    for prs_idx, prs in enumerate(atm.p):
        if prs < trpp_prs:
            atm.tmp[prs_idx] = trpp_tmp

    # Staggered nodes
    for prsl_idx, prls in enumerate(atm.pl):
        if prls < trpp_prs:
            atm.tmpl[prsl_idx] = trpp_tmp

    # Set mixing ratios to same as tropopause
    for idx in reversed(range(0, trpp_idx)):
    
        atm.cp[idx] = 0.

        # Volatile abundances
        for vol in atm.vol_list.keys():
        
            # Saturation vapor pressure
            p_vol_sat     = ga.p_sat(vol, atm.tmp[idx])

            # If still condensible
            if atm.p[idx] > p_vol_sat:

                cond_diff            = (atm.x_cond[vol][idx] - atm.x_cond[vol][trpp_idx])

                atm.xc[idx]          -= cond_diff
                atm.xv[idx]          += cond_diff

                atm.x_cond[vol][idx] = atm.x_cond[vol][trpp_idx]
                atm.x_gas[vol][idx]  = atm.x_gas[vol][trpp_idx]
                atm.p_vol[vol][idx]  = atm.x_gas[vol][idx] * atm.p[idx]

            # If not anymore
            else:
                atm.xc[idx]          -= atm.x_cond[vol][idx]
                atm.x_gas[vol][idx]  = atm.x_gas[vol][trpp_idx]
                atm.xd[idx]          += atm.x_gas[vol][idx] + atm.x_cond[vol][idx]
                atm.xv[idx]          -= atm.x_gas[vol][idx]

                atm.x_cond[vol][idx] = 0.
                
                atm.p_vol[vol][idx]  = atm.x_gas[vol][trpp_idx] * atm.p[idx]
                

            # Renormalize cp w/ molar concentration
            atm.cp[idx]   += (atm.x_gas[vol][idx] + atm.x_cond[vol][idx]) * ga.cpv(vol, atm.tmp[idx]) / (atm.xd[idx] + atm.xv[idx] + atm.xc[idx]) # w/ cond

    return atm

def InterpolateStellarLuminosity(star_mass, star_age, mean_distance):

    luminosity_df           = pd.read_csv("luminosity_tracks/Lum_m"+str(star_mass)+".txt")
    star_age                = star_age/1e+6         # Myr
    ages                    = luminosity_df["age"]*1e+3 # Myr
    luminosities            = luminosity_df["lum"]      # L_sol

    # Interpolate luminosity for current time
    interpolate_luminosity  = interpolate.interp1d(ages, luminosities)
    interpolated_luminosity = interpolate_luminosity([star_age])*L_sun

    stellar_toa_heating     = interpolated_luminosity / ( 4. * np.pi * (mean_distance*AU)**2. )
    
    return stellar_toa_heating[0], interpolated_luminosity/L_sun

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    ##### Settings

    # Planet age and orbit
    star_age      = 4567e+6              # yr
    mean_distance = 1.0                 # au

    # Surface pressure & temperature
    P_surf        = 10e+5              # Pa
    T_surf        = 300.               # K

    # Volatile molar concentrations: ! must sum to one !
    vol_list = { 
                  "H2O" : .9999, 
                  "CO2" : .0,
                  "H2"  : .0, 
                  "N2"  : .0,  
                  "CH4" : .0, 
                  "O2"  : .0, 
                  "CO"  : .0, 
                  "He"  : .0,
                  "NH3" : .0, 
                }

    # Stellar heating on/off
    stellar_heating = True

    # Compare to the dry adiabat w/ time stepping
    cp_dry = False

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, vol_list)

    # Compute stellar heating
    toa_heating, star_luminosity = InterpolateStellarLuminosity(1.0, star_age, mean_distance)

    # Set stellar heating on or off
    if stellar_heating == False: toa_heating = 0.

    # Compute heat flux
    atm = RadConvEqm("./output", star_age, atm, toa_heating, [], [], standalone=True, cp_dry=False)

    # Plot abundances w/ TP structure
    ga.plot_adiabats(atm)
