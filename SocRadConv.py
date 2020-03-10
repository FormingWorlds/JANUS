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

# Parameters to run SocRadConv stand-alone
rad_steps = 100

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

# Calculate dew point temperature
def Tdew_H2O(p):
    return Tref/(1-(Tref*R/L)*math.log(p/pref))

# Moist adjustment switch
Moist_Adjustment = True

# Number of dry adjustment steps
Nsteps_dry = 5

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

    # Set up pressure array (a global)
    pstart          = atm.ps#*.995
    rat             = (atm.ptop/pstart)**(1./atm.nlev)
    logLevels       = [pstart*rat**i for i in range(atm.nlev+1)]
    logLevels.reverse()
    levels          = [atm.ptop + i*(pstart-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl          = np.array(logLevels)
    atm.p           = (atm.pl[1:] + atm.pl[:-1]) / 2

    # Initialize on adiabat
    atm.temp        = atm.ts*(atm.p/atm.p[-1])**atm.Rcp

    # Set initial stratosphere guess to isothermal (closer to actual solution)
    atm.temp        = np.where(atm.temp<atm.ts/4.,atm.ts/4.,atm.temp) 

    # Calculate individual moist adiabats
    Moist_adiabat_H2O   = [ Tdew_H2O(pp) for pp in atm.p ]

    # TdewH2O = [ ga.Tdew( 'H2O', pressure ) for pressure in atm.p ]
    # TdewCO2 = [ ga.Tdew( 'CO2', pressure ) for pressure in atm.p ]
    # TdewCH4 = [ ga.Tdew( 'CH4', pressure ) for pressure in atm.p ]
    # TdewCO  = [ ga.Tdew( 'CO',  pressure ) for pressure in atm.p ]
    # TdewN2  = [ ga.Tdew( 'N2',  pressure ) for pressure in atm.p ]
    # TdewO2  = [ ga.Tdew( 'O2',  pressure ) for pressure in atm.p ]
    # TdewH2  = [ ga.Tdew( 'H2',  pressure ) for pressure in atm.p ]
    # TdewHe  = [ ga.Tdew( 'He',  pressure ) for pressure in atm.p ]
    # TdewNH3 = [ ga.Tdew( 'NH3', pressure ) for pressure in atm.p ]
    
    # Feed mixing ratios
    if standalone == True:
        atm.mixing_ratios[0] = atm_chemistry["H2O"]    # H2O
        atm.mixing_ratios[1] = atm_chemistry["CO2"]    # CO2
        atm.mixing_ratios[2] = atm_chemistry["H2"]     # H2
        atm.mixing_ratios[3] = atm_chemistry["CH4"]    # CH4
        atm.mixing_ratios[4] = atm_chemistry["CO"]     # CO
        atm.mixing_ratios[5] = atm_chemistry["N2"]     # N2
        atm.mixing_ratios[6] = atm_chemistry["O2"]     # O2

        use_vulcan = 0
    else:
        # Varying mixing ratios with height
        if SPIDER_options["use_vulcan"] == 1 or SPIDER_options["use_vulcan"] == 2:
            atm_chemistry = atm_chemistry.reindex(index=atm_chemistry.index[::-1])
            atm.mixing_ratios[0] = atm_chemistry["H2O"]    # H2O
            atm.mixing_ratios[1] = atm_chemistry["CO2"]    # CO2
            atm.mixing_ratios[2] = atm_chemistry["H2"]     # H2
            atm.mixing_ratios[3] = atm_chemistry["CH4"]    # CH4
            atm.mixing_ratios[4] = atm_chemistry["CO"]     # CO
            atm.mixing_ratios[5] = atm_chemistry["N2"]     # N2
            atm.mixing_ratios[6] = atm_chemistry["O2"]     # O2
        # Constant mixing ratios with SOCRATES
        else: 
            atm.mixing_ratios[0] = runtime_helpfile.iloc[-1]["H2O_mr"]    # H2O
            atm.mixing_ratios[1] = runtime_helpfile.iloc[-1]["CO2_mr"]    # CO2
            atm.mixing_ratios[2] = runtime_helpfile.iloc[-1]["H2_mr"]     # H2
            atm.mixing_ratios[3] = runtime_helpfile.iloc[-1]["CH4_mr"]    # CH4
            atm.mixing_ratios[4] = runtime_helpfile.iloc[-1]["CO_mr"]     # CO
            atm.mixing_ratios[5] = runtime_helpfile.iloc[-1]["N2_mr"]     # N2
            atm.mixing_ratios[6] = runtime_helpfile.iloc[-1]["O2_mr"]     # O2

        use_vulcan = SPIDER_options["use_vulcan"]

    # Initialise previous OLR and TOA heating to zero
    PrevOLR     = 0.
    PrevMaxHeat = 0.
    PrevTemp    = 0.*atm.temp[:]

    # Initialization complete
    # Now do the time stepping
    # matplotlib.rc('axes',edgecolor='k')
    for i in range(0,rad_steps):

        atm, moist_wo_cond, moist_w_cond = steps(atm, stellar_toa_heating, atm_chemistry, use_vulcan)

        if i % 1 == 0:

            sns.set_style("ticks")
            sns.despine()
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

            ax1.semilogy(atm.temp,atm.p*1e-5, ls="-", label=r'Dry adiabat')

            if Moist_Adjustment == True:
                # ax1.semilogy(moist_wo_cond, atm.p*1e-5, color="green", ls="-", label=r'Moist adiabat, no cond.')

                # ax1.semilogy(moist_w_cond, atm.p*1e-5, color="green", ls="--", label=r'Moist adiabat')
                
                # ax1.semilogy(Moist_adiabat_H2O, atm.p*1e-5, color="red", ls="-.", label=r'H$_2$O moist adiabat')

                ax1.semilogy(TdewH2O, atm.p*1e-5, label=r'H$_2$O', lw=0.8, ls=":")
                ax1.semilogy(TdewCO2, atm.p*1e-5, label=r'CO$_2$', lw=0.8, ls=":")
                ax1.semilogy(TdewH2,  atm.p*1e-5, label=r'H$_2$',  lw=0.8, ls=":")
                ax1.semilogy(TdewCH4, atm.p*1e-5, label=r'CH$_4$', lw=0.8, ls=":")
                ax1.semilogy(TdewCO,  atm.p*1e-5, label=r'CO',     lw=0.8, ls=":")
                ax1.semilogy(TdewN2,  atm.p*1e-5, label=r'N$_2$',  lw=0.8, ls=":")
                ax1.semilogy(TdewO2,  atm.p*1e-5, label=r'O$_2$',  lw=0.8, ls=":")

            ax1.invert_yaxis()
            ax1.set_xlabel('Temperature (K)')
            ax1.set_ylabel('Pressure (bar)')
            ax1.set_xlim([0,np.max(atm.temp)])
            ax1.set_ylim([np.max(atm.p*1e-5),np.min(atm.p*1e-5)])
            ax1.set_xticks([0,0.2*np.max(atm.temp),0.4*np.max(atm.temp),0.6*np.max(atm.temp),0.8*np.max(atm.temp),np.max(atm.temp)])
            ax1.legend()
            
            ax2.plot(atm.band_centres,atm.LW_spectral_flux_up[:,0]/atm.band_widths)
            ax2.plot(atm.band_centres,surf_Planck_nu(atm)/atm.band_widths,color="gray",ls='--',label='Black body ('+str(atm.ts)+" K)")
            ax2.set_xlim([np.min(atm.band_centres),np.max(atm.band_centres)])
            ax2.set_ylabel('Spectral flux density (Jy?)')
            ax2.set_xlabel('Wavenumber (1/cm)')
            ax2.legend()

            plt.savefig(output_dir+'/TP_profile_'+str(round(time_current))+'.pdf', bbox_inches="tight")
            plt.close(fig)
            print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)))

        if i % 10 == 0:
            print("Iteration", i, end =", ")
            print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)), ", dt =", atm.dt)

        # Reduce timestep if heating is not converging
        if abs(np.max(atm.temp-PrevTemp[:])) < 0.05 or abs(atm.temp[0]-atm.temp[1]) > 3.0:
            atm.dt  = atm.dt*0.99
            # print("Not converging -> reduce timestep to dt =", atm.dt)

        # Sensitivity break condition
        if (abs(atm.LW_flux_up[0]-PrevOLR) < (0.1*(5.67e-8*atm.ts**4)**0.5)) and i > 5 :
           print("Break -> deltaOLR =", abs(atm.LW_flux_up[0]-PrevOLR), ", deltaT =", abs(np.max(atm.temp-PrevTemp[:])))
           break    # break here

        PrevOLR     = atm.LW_flux_up[0]
        PrevMaxHeat = abs(np.max(atm.total_heating))
        PrevTemp[:] = atm.temp[:]

    # Write TP and spectral flux profiles for later plotting
    out_a = np.column_stack( ( atm.temp, atm.p*1e-5 ) ) # K, Pa->bar
    np.savetxt( output_dir+str(int(time_current))+"_atm_TP_profile.dat", out_a )
    out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths ) )
    np.savetxt( output_dir+str(int(time_current))+"_atm_spectral_flux.dat", out_a )

    return atm.LW_flux_up[-1], atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths
  
# Dry convective adjustment routine
def dryAdj(atm):

    T = atm.temp
    p = atm.p
    
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
            atm.temp[i] = T1
            atm.temp[i+1] = T2
    
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
            atm.temp[i]   = T1
            atm.temp[i+1] = T2       

# # Moist adjustment routine
# def moistAdj(atm, atm_chemistry, use_vulcan):

#     moist_adiabat = ga.solve_general_adiabat(atm, atm_chemistry, use_vulcan)
        
#     # moist_adiabat = ga.General_moist_adiabat(atm.ps,T,xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3)

#     # print(Tdew)

#     # # Downward pass
#     # for i in range(len(T)-1):
#     #     if T[i] < Tdew(p[i]):
#     #         T[i] = Tdew(p[i])   # Temperature stays the same during phase change
#     # # Upward pass
#     # for i in range(len(T)-2,-1,-1): 
#     #     if T[i] < Tdew(p[i]):
#     #         T[i] = Tdew(p[i])

#     return moist_adiabat

# Time integration for n steps
def steps(atm, stellar_toa_heating, atm_chemistry, use_vulcan):
    atm     = SocRadModel.radCompSoc(atm, stellar_toa_heating)
    dT      = atm.total_heating*atm.dt
    
    # Limit the temperature change per step
    dT      = np.where(dT>5.,5.,dT)
    dT      = np.where(dT<-5.,-5.,dT)
    
    # Midpoint method time stepping
    # changed call to r.  Also modified to hold Tg fixed
    atm     = SocRadModel.radCompSoc(atm, stellar_toa_heating)
    dT      = atm.total_heating*atm.dt
    
    # Limit the temperature change per step
    dT      = np.where(dT>5.,5.,dT)
    dT      = np.where(dT<-5.,-5.,dT)
    atm.temp += dT
    dTmax   = max(abs(dT)) #To keep track of convergence

    # Do the surface balance
    kturb   = .1
    atm.temp[-1] += -atm.dt*kturb*(atm.temp[-1] - atm.ts)
    
    # Adiabatic adjustment
    for iadj in range(Nsteps_dry):
        
        # Dry adjustment step
        dryAdj(atm)

    # Moist adjustment step
    if Moist_Adjustment == True:
        atm = ga.general_adiabat(copy.deepcopy(atm))
        # print(moist_wo_cond, moist_w_cond)
    else:
        moist_wo_cond, moist_w_cond = np.zeros(len(atm.temp))

    Tad = atm.temp[-1]*(atm.p/atm.p[-1])**atm.Rcp
    
    return atm, moist_wo_cond, moist_w_cond

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

##### Stand-alone initial conditions

# Planet age and orbit
time_current  = 1e+7                # yr
time_offset   = 1e+7                # yr
mean_distance = 1.0                 # au

# Surface pressure & temperature
P_surf                  = 1e+3      # Pa
T_surf                  = 280.      # K

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

# Define the atmosphere object from this input

# Create atmosphere object
atm                     = atmos()

# Surface temperature of planet
atm.ts                  = T_surf     # K
atm.tmp[0]              = atm.ts  

# Surface & top pressure
atm.ps                  = P_surf     # Pa
atm.p[0]                = atm.ps     # Pa
atm.ptop                = np.amin([atm.ps*1e-10, 1e-5])   # Pa

# Initialize level-dependent quantities
atm.vol_list            = vol_list   # List of all species and initial concentrations

# Compute stellar heating
toa_heating, star_luminosity = InterpolateStellarLuminosity(1.0, time_current, time_offset, mean_distance)

# Construct radiative-convective profile + heat flux
LW_flux_up, band_centres, LW_spectral_flux_up_per_band_widths = RadConvEqm("./output", time_current, atm, toa_heating, vol_list, [], [], standalone=True)
