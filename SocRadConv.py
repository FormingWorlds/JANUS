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



# Thermodynamic constants for moist adjustment
R         = phys.water.R                  # J/(kg*K) specific gas constant of water vapor
Rcp       = phys.water.Rcp                # cp in J/(kg*K) specific heat constant of water vapor
L         = phys.water.L_vaporization     # J/kg, latent heat of condensation of water vapor at 300K
# Saturation vapor pressure, args=T,T0,e0,MolecularWeight,LatentHeat
esat      = phys.satvps_function(phys.water) 
Tref      = 2000.                         # Reference temperature
pref      = esat(Tref)                    # Reference pressure

# Other constants
L_sun     = 3.828e+26                     # W, IAU definition
AU        = 1.495978707e+11               # m

R_universal = 8.31446261815324            # Universal gas constant, J.K-1.mol-1

# # Calculate dew point temperature
# def Tdew_H2O(p):
#     return Tref/(1-(Tref*R/L)*math.log(p/pref))

# Number of radiation steps
rad_steps  = 1

# Number of convective adjustment steps
conv_steps = 0

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

    # # Instantiate the radiation model
    # atm = atmos()

    # if standalone == True:
    #     atm.ps      = runtime_helpfile[0]*1e5   # bar->Pa
    #     atm.ts      = runtime_helpfile[1]       # K   
    # else:
    #     atm.ps      = runtime_helpfile.iloc[-1]["P_surf"]*1e5 # bar->Pa
    #     atm.ts      = runtime_helpfile.iloc[-1]["T_surf"]

    # # Avoid math error for moist adiabat
    # atm.ptop        = atm.ps*1e-5   

    # # Set up pressure array (a global)
    # pstart          = atm.ps#*.995
    # rat             = (atm.ptop/pstart)**(1./atm.nlev)
    # logLevels       = [pstart*rat**i for i in range(atm.nlev+1)]
    # logLevels.reverse()
    # levels          = [atm.ptop + i*(pstart-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    # plevels2          = np.array(logLevels)
    # plevels           = (plevels2[1:] + plevels2[:-1]) / 2

    # print(atm.ps)
    # print(plevels2)
    # print(plevels)

    # # Initialize on general adiabat
    # atm = ga.general_adiabat(atm)

    # print(len(atm.p), len(atm.pl), len(atm.tmp), len(atm.tmpl))
    # for idx, val in enumerate(atm.p):
    #     print(atm.p[idx], atm.pl[idx], atm.tmp[idx], atm.tmpl[idx])
    # # Initialize on adiabat
    # atm.tmp        = atm.ts*(atm.p/atm.p[-1])**atm.Rcp

    # # Set initial stratosphere guess to isothermal (closer to actual solution)
    # atm.tmp        = np.where(atm.temp<atm.ts/4.,atm.ts/4.,atm.temp) 

    # # Calculate individual moist adiabats
    # Moist_adiabat_H2O   = [ Tdew_H2O(pp) for pp in atm.p ]

    # TdewH2O = [ ga.Tdew( 'H2O', pressure ) for pressure in atm.p ]
    # TdewCO2 = [ ga.Tdew( 'CO2', pressure ) for pressure in atm.p ]
    # TdewCH4 = [ ga.Tdew( 'CH4', pressure ) for pressure in atm.p ]
    # TdewCO  = [ ga.Tdew( 'CO',  pressure ) for pressure in atm.p ]
    # TdewN2  = [ ga.Tdew( 'N2',  pressure ) for pressure in atm.p ]
    # TdewO2  = [ ga.Tdew( 'O2',  pressure ) for pressure in atm.p ]
    # TdewH2  = [ ga.Tdew( 'H2',  pressure ) for pressure in atm.p ]
    # TdewHe  = [ ga.Tdew( 'He',  pressure ) for pressure in atm.p ]
    # TdewNH3 = [ ga.Tdew( 'NH3', pressure ) for pressure in atm.p ]
    
    # # Feed mixing ratios
    # if standalone == True:
    #     atm.mixing_ratios[0] = atm_chemistry["H2O"]    # H2O
    #     atm.mixing_ratios[1] = atm_chemistry["CO2"]    # CO2
    #     atm.mixing_ratios[2] = atm_chemistry["H2"]     # H2
    #     atm.mixing_ratios[3] = atm_chemistry["CH4"]    # CH4
    #     atm.mixing_ratios[4] = atm_chemistry["CO"]     # CO
    #     atm.mixing_ratios[5] = atm_chemistry["N2"]     # N2
    #     atm.mixing_ratios[6] = atm_chemistry["O2"]     # O2

    #     use_vulcan = 0
    # else:
    #     # Varying mixing ratios with height
    #     if SPIDER_options["use_vulcan"] == 1 or SPIDER_options["use_vulcan"] == 2:
    #         atm_chemistry = atm_chemistry.reindex(index=atm_chemistry.index[::-1])
    #         atm.mixing_ratios[0] = atm_chemistry["H2O"]    # H2O
    #         atm.mixing_ratios[1] = atm_chemistry["CO2"]    # CO2
    #         atm.mixing_ratios[2] = atm_chemistry["H2"]     # H2
    #         atm.mixing_ratios[3] = atm_chemistry["CH4"]    # CH4
    #         atm.mixing_ratios[4] = atm_chemistry["CO"]     # CO
    #         atm.mixing_ratios[5] = atm_chemistry["N2"]     # N2
    #         atm.mixing_ratios[6] = atm_chemistry["O2"]     # O2
    #     # Constant mixing ratios with SOCRATES
    #     else: 
    #         atm.mixing_ratios[0] = runtime_helpfile.iloc[-1]["H2O_mr"]    # H2O
    #         atm.mixing_ratios[1] = runtime_helpfile.iloc[-1]["CO2_mr"]    # CO2
    #         atm.mixing_ratios[2] = runtime_helpfile.iloc[-1]["H2_mr"]     # H2
    #         atm.mixing_ratios[3] = runtime_helpfile.iloc[-1]["CH4_mr"]    # CH4
    #         atm.mixing_ratios[4] = runtime_helpfile.iloc[-1]["CO_mr"]     # CO
    #         atm.mixing_ratios[5] = runtime_helpfile.iloc[-1]["N2_mr"]     # N2
    #         atm.mixing_ratios[6] = runtime_helpfile.iloc[-1]["O2_mr"]     # O2

    #     use_vulcan = SPIDER_options["use_vulcan"]

    # # Initialise previous OLR and TOA heating to zero
    # PrevOLR     = 0.
    # PrevMaxHeat = 0.
    # PrevTemp    = 0.*atm.tmp

    # Initialization complete
    # Now do the time stepping
    # matplotlib.rc('axes',edgecolor='k')
    # for i in range(0, rad_steps):

    atm_dry, atm_moist = radiation_timestepping(atm, toa_heating, rad_steps)

    # if i % 1 == 0:

    #     sns.set_style("ticks")
    #     sns.despine()
    #     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

    #     ax1.semilogy(atm.temp,atm.p*1e-5, ls="-", label=r'Dry adiabat')

    #     ax1.invert_yaxis()
    #     ax1.set_xlabel('Temperature (K)')
    #     ax1.set_ylabel('Pressure (bar)')
    #     ax1.set_xlim([0,np.max(atm.temp)])
    #     ax1.set_ylim([np.max(atm.p*1e-5),np.min(atm.p*1e-5)])
    #     ax1.set_xticks([0,0.2*np.max(atm.temp),0.4*np.max(atm.temp),0.6*np.max(atm.temp),0.8*np.max(atm.temp),np.max(atm.temp)])
    #     ax1.legend()
        
    #     ax2.plot(atm.band_centres,atm.LW_spectral_flux_up[:,0]/atm.band_widths)
    #     ax2.plot(atm.band_centres,surf_Planck_nu(atm)/atm.band_widths,color="gray",ls='--',label='Black body ('+str(atm.ts)+" K)")
    #     ax2.set_xlim([np.min(atm.band_centres),np.max(atm.band_centres)])
    #     ax2.set_ylabel('Spectral flux density (Jy?)')
    #     ax2.set_xlabel('Wavenumber (1/cm)')
    #     ax2.legend()

    #     plt.savefig(output_dir+'/TP_profile_'+str(round(time_current))+'.pdf', bbox_inches="tight")
    #     plt.close(fig)
    #     print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)))

        # PrevOLR     = atm.LW_flux_up[0]
        # PrevMaxHeat = abs(np.max(atm.total_heating))
        # PrevTemp[:] = atm.temp[:]

    # # Plot results
    # plot_heat_balance(atm)

    plot_heat_balance(atm_dry, atm_moist)

    # Write TP and spectral flux profiles for later plotting
    out_a = np.column_stack( ( atm.temp, atm.p*1e-5 ) ) # K, Pa->bar
    np.savetxt( output_dir+str(int(time_current))+"_atm_TP_profile.dat", out_a )
    out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths ) )
    np.savetxt( output_dir+str(int(time_current))+"_atm_spectral_flux.dat", out_a )

    return atm.LW_flux_up[-1], atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths


# Dry adiabat profile
def dry_adiabat_atm(atm):

    # # Define pressure levels
    # atm     = set_pressure_array(atm)

    # # Define pressure levels
    # atm.     = set_pressure_array(atm)

    # Calculate cp from molar concentrations
    cp_mix = 0.
    for vol in atm.vol_list.keys():
        cp_mix += atm.vol_list[vol] * ga.cpv(vol)

    # Calculate dry adiabat slope
    atm.Rcp = R_universal / cp_mix

    # Calculate dry adiabat temperature profile
    for idx, prs in enumerate(atm.p):
        atm.tmp[idx] = atm.ts * ( atm.p[idx] / atm.ps ) ** ( atm.Rcp )

    atm.tmpl = np.interp(atm.pl, np.flip(atm.p), np.flip(atm.tmp))

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
            atm.tmp[i] = T1
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

        ax1.invert_yaxis()
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Pressure (Pa)')

        # ax1.set_xlim([0,np.max(atm.temp)])
        # ax1.set_ylim([np.max(atm.p*1e-5),np.min(atm.p*1e-5)])

        # ax1.set_xticks([0,0.2*np.max(atm.temp),0.4*np.max(atm.temp),0.6*np.max(atm.temp),0.8*np.max(atm.temp),np.max(atm.temp)])

        ax1.legend()
        
        ax2.plot(atm_moist.band_centres,surf_Planck_nu(atm_moist)/atm_moist.band_widths,color="gray",ls='--',label='Black body ('+str(atm_moist.ts)+" K)")

        # ax2.plot(atm_dry.band_centres,atm_dry.LW_spectral_flux_up[:,0]/atm_dry.band_widths)
        # ax2.plot(atm_moist.band_centres,atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths)
        
        # ax2.set_xlim([np.min(atm.band_centres),np.max(atm.band_centres)])

        ax2.set_ylabel('Spectral flux density (Jy?)')
        ax2.set_xlabel('Wavenumber (1/cm)')
        ax2.legend()

        plt.savefig("./output"+'/TP_profile_'+str(round(time_current))+'.pdf', bbox_inches="tight")
        plt.close(fig)
        # print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)))

# Time integration for n steps
def radiation_timestepping(atm, toa_heating, rad_steps):

    # Initialise previous OLR and TOA heating to zero
    PrevOLR_dry         = 0.
    PrevMaxHeat_dry     = 0.
    PrevTemp_dry        = atm.tmp * 0.
    PrevOLR_moist       = 0.
    PrevMaxHeat_moist   = 0.
    PrevTemp_moist      = atm.tmp * 0.

    # # Create deep copies that are distinct from old dict
    # atm_dry             = copy.deepcopy(atm)
    # atm_moist           = copy.deepcopy(atm)

    # Build moist adiabat structure
    atm_moist           = ga.general_adiabat(atm)
    atm_dry             = dry_adiabat_atm(copy.deepcopy(atm_moist))
    
    plot_heat_balance(atm_dry, atm_moist)

    # Time stepping
    for i in range(0, rad_steps):

        ### Dry calculation

        print(len(atm_dry.tmp), atm_dry.tmp)

        # Compute radiation, midpoint method time stepping
        atm_dry         = SocRadModel.radCompSoc(atm_dry, toa_heating)
        dT_dry          = atm_dry.total_heating * atm_dry.dt

        # Limit the temperature change per step
        dT_dry          = np.where(dT_dry > 5., 5., dT_dry)
        dT_dry          = np.where(dT_dry < -5., -5., dT_dry)

        print(len(atm_dry.tmp), atm_dry.tmp)
        print(len(dT_dry), dT_dry)

        atm_dry.tmp     += dT_dry

        # To keep track of convergence
        dTmax_dry       = max(abs(dT_dry)) 

        # # Do the surface balance
        # kturb       = .1
        # atm.tmp[-1] += -atm.dt * kturb * (atm.tmp[-1] - atm.ts)
        
        # Dry adiabatic adjustment
        for iadj in range(conv_steps):
            atm_dry     = DryAdj(atm_dry)

        # plot_heat_balance(atm_dry, atm_moist)

        ### Moist calculation

        # Compute radiation, midpoint method time stepping
        atm_moist       = SocRadModel.radCompSoc(atm_moist, toa_heating)
        dT_moist        = atm_moist.total_heating * atm_moist.dt

        # Limit the temperature change per step
        dT_moist        = np.where(dT_moist > 5., 5., dT_moist)
        dT_moist        = np.where(dT_moist < -5., -5., dT_moist)

        # To keep track of convergence
        dTmax_moist     = max(abs(dT_moist)) 

        # Moist single-step adjustment
        atm_moist       = MoistAdj(atm_moist, dT_moist)

        # Tad = atm.tmp[-1]*(atm.p/atm.p[-1])**atm.Rcp

        # Inform during runtime
        if i % 10 == 0:
            print("Iteration ", i, end=" (dry, moist): ")
            print("OLR: " + str(PrevOLR_dry) + ", " + str(PrevOLR_moist) + " W/m^2,", "dT_max = " + str(dTmax_dry) + ", " + str(dTmax_moist) + " K")

        # Reduce timestep if heating is not converging
        dTglobal = abs(np.max(atm_dry.tmp-PrevTemp_dry[:]))
        dTtop    = abs(atm_dry.tmp[0]-atm_dry.tmp[1])
        if dTglobal < 0.05 or dTtop > 3.0:
            atm_dry.dt  = atm_dry.dt*0.99
            print("Dry adiabat structure not converging -> dt_new =", atm_dry.dt)
        dTglobal = abs(np.max(atm_moist.tmp-PrevTemp_moist[:]))
        dTtop    = abs(atm_moist.tmp[0]-atm_moist.tmp[1])
        if dTglobal < 0.05 or dTtop > 3.0:
            atm_moist.dt  = atm_moist.dt*0.99
            print("Moist adiabat structure not converging -> dt_new =", atm_moist.dt)

        # Sensitivity break condition
        if (abs(atm_dry.LW_flux_up[0]-PrevOLR_dry) < (0.1*(5.67e-8*atm_dry.ts**4)**0.5)) and i > 5 :
           print("Break -> deltaOLR =", abs(atm_dry.LW_flux_up[0]-PrevOLR_dry), ", deltaT =", abs(np.max(atm_dry.tmp-PrevTemp_dry[:])))
           break    # break here
        if (abs(atm_moist.LW_flux_up[0]-PrevOLR_moist) < (0.1*(5.67e-8*atm_moist.ts**4)**0.5)) and i > 5 :
           print("Break -> deltaOLR =", abs(atm_moist.LW_flux_up[0]-PrevOLR_moist), ", deltaT =", abs(np.max(atm_moist.tmp-PrevTemp_moist[:])))
           break    # break here

        PrevOLR_dry       = atm_dry.LW_flux_up[0]
        PrevMaxHeat_dry   = abs(np.max(atm_dry.total_heating))
        PrevTemp_dry[:]   = atm_dry.tmp[:]
        PrevOLR_moist     = atm_moist.LW_flux_up[0]
        PrevMaxHeat_moist = abs(np.max(atm_moist.total_heating))
        PrevTemp_moist[:] = atm_moist.tmp[:]
    
    return atm_dry, atm_moist

# Define pressure levels for dry adjustment
def set_pressure_array(atm):
   
    # MHD pressure levels
    # rat       = (atm.ptop/atm.ps)**(1./atm.nlev)
    # logLevels = [atm.ps*rat**i for i in range(atm.nlev+1)]
    # levels    = [atm.ptop + i*(atm.ps-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev)]

    # Logspace pressure levels
    atm.p     = np.flip(np.logspace(np.log10(atm.ptop), np.log10(atm.ps), atm.nlev))

    # Interpolate staggered nodes
    for idx, prs in enumerate(atm.p):
        print(idx, prs)
        if idx == 0:
            atm.pl[idx] = atm.p[idx]
        else:
            atm.pl[idx] = (atm.p[idx-1] + atm.p[idx]) / 2.
        if idx == len(atm.p):
            atm.pl[-1] = atm.p[idx] - ((atm.p[-2]-atm.p[-1])/2.)

    # atm.pl    = (atm.p[1:] + atm.p[:-1]) / 2

    return atm

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
P_surf        = 1e+5                # Pa
T_surf        = 800.                # K

# Volatile molar concentrations: ! must sum to one !
vol_list = { 
              "H2O" : .5, 
              "CO2" : .5,
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
