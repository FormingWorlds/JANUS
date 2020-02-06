'''
MDH 28/01/19
Socrates radiative-convective model
'''

import numpy as np
import math,phys
import General_adiabat as ga # Moist adiabat with multiple condensibles
import matplotlib.pyplot as plt
import matplotlib
import SocRadModel
from atmosphere_column import atmos
import pandas as pd

# Parameters to run SocRadConv stand-alone
rad_steps = 100

# Thermodynamic constants for moist adjustment
R    = phys.water.R                  # J/(kg*K) specific gas constant of water vapor
Rcp  = phys.water.Rcp                # cp in J/(kg*K) specific heat constant of water vapor
L    = phys.water.L_vaporization     # J/kg, latent heat of condensation of water vapor at 300K
# Saturation vapor pressure, args=T,T0,e0,MolecularWeight,LatentHeat
esat = phys.satvps_function(phys.water) 
Tref = 3000.                         # Reference temperature
pref = esat(Tref)                    # Reference pressure

# Calculate dew point temperature
def Tdew(p):
    return Tref/(1-(Tref*R/L)*math.log(p/pref))

# Moist adjustment switch
Moist_Adjustment = False

def surf_Planck_nu(atm):
    h  = 6.63e-34
    c  = 3.0e8
    kb = 1.38e-23
    B  = np.zeros(len(atm.band_centres))
    c1 = 1.191042e-5
    c2 = 1.4387752
    for i in range(len(atm.band_centres)):
        nu = atm.band_centres[i]
        B[i] = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))

    B = B * atm.band_widths/1000.0
    return B

def RadConvEqm(output_dir, time_current, runtime_helpfile, stellar_toa_heating, atm_chemistry, loop_counter, SPIDER_options):

    # Instantiate the radiation model
    atm = atmos()

    # Set up pressure array (a global)
    atm.ps      = runtime_helpfile.iloc[-1]["P_surf"]*1e5 # bar->Pa
    pstart      = atm.ps#*.995
    rat         = (atm.ptop/pstart)**(1./atm.nlev)
    logLevels   = [pstart*rat**i for i in range(atm.nlev+1)]
    logLevels.reverse()
    levels      = [atm.ptop + i*(pstart-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl      = np.array(logLevels)
    atm.p       = (atm.pl[1:] + atm.pl[:-1]) / 2


    # Now do the calculation
    atm.ts          = runtime_helpfile.iloc[-1]["T_surf"]   
    atm.Rcp         = 2./7.
    atm.temp        = atm.ts*(atm.p/atm.p[-1])**atm.Rcp               # Initialize on adiabat
    atm.temp        = np.where(atm.temp<atm.ts/4.,atm.ts/4.,atm.temp) # CHECK
    # atm.temp  = np.where(atm.temp<400.,400.,atm.temp)
    atm.n_species   = 7
    Moist_adiabat   = [ Tdew(pp) for pp in atm.p ]
    
    # Feed mixing ratios
    if SPIDER_options["use_vulcan"] == 1 or SPIDER_options["use_vulcan"] == 2:
        atm_chemistry = atm_chemistry.reindex(index=atm_chemistry.index[::-1])
        atm.mixing_ratios[0] = atm_chemistry["H2O"]    # H2O
        atm.mixing_ratios[1] = atm_chemistry["CO2"]    # CO2
        atm.mixing_ratios[2] = atm_chemistry["H2"]     # H2
        atm.mixing_ratios[3] = atm_chemistry["CH4"]    # CH4
        atm.mixing_ratios[4] = atm_chemistry["CO"]     # CO
        atm.mixing_ratios[5] = atm_chemistry["N2"]     # N2
        atm.mixing_ratios[6] = atm_chemistry["O2"]     # O2
    else: # Constant mixing ratios from SPIDER
        atm.mixing_ratios[0] = runtime_helpfile.iloc[-1]["H2O_mr"]    # H2O
        atm.mixing_ratios[1] = runtime_helpfile.iloc[-1]["CO2_mr"]    # CO2
        atm.mixing_ratios[2] = runtime_helpfile.iloc[-1]["H2_mr"]     # H2
        atm.mixing_ratios[3] = runtime_helpfile.iloc[-1]["CH4_mr"]    # CH4
        atm.mixing_ratios[4] = runtime_helpfile.iloc[-1]["CO_mr"]     # CO
        atm.mixing_ratios[5] = runtime_helpfile.iloc[-1]["N2_mr"]     # N2
        atm.mixing_ratios[6] = runtime_helpfile.iloc[-1]["O2_mr"]     # O2

    # Initialise previous OLR and TOA heating to zero
    PrevOLR     = 0.
    PrevMaxHeat = 0.
    PrevTemp    = 0.*atm.temp[:]

    # Initialization complete
    # Now do the time stepping
    # matplotlib.rc('axes',edgecolor='k')
    for i in range(0,rad_steps):

        atm = steps(atm, stellar_toa_heating)

        #hack!
        # atm.temp[0] = atm.temp[1]

        if i % 1 == 0:
            if 1==2:
                plt.figure(figsize=(7,4))
                plt.semilogy(atm.temp,atm.p)
                plt.gca().invert_yaxis()
                plt.ylabel('Pressure [mb]')
                plt.xlabel('Temperature [K]')
                plt.gca().xaxis.label.set_color('white')
                plt.tick_params(axis='x', colors='white')
                plt.gca().yaxis.label.set_color('white')
                plt.tick_params(axis='y', colors='white')
                # plt.show()
            #print("OLR " + str(atm.LW_flux_up[-1]))
            #print("OLR change " + str(atm.LW_flux_up[-1]-PrevOLR))
            # print("Max heating " + str(np.max(atm.total_heating)))
            #print("Max dT " + str(abs(np.max(atm.temp-PrevTemp[:]))))
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,7))
            ax1.semilogy(atm.temp,atm.p)
            ax1.semilogy(Moist_adiabat,atm.p,'r',label='Moist adiabat')
            ax1.invert_yaxis()
            ax1.set_xlabel('Temperature [K]')
            ax1.set_ylabel('Pressure [mb]')
            ax2.plot(atm.band_centres,atm.LW_spectral_flux_up[:,0]/atm.band_widths,'k')
#            ax2.plot(atm.band_centres,surf_Planck_nu(atm)/atm.band_widths,'k--')
            ax2.set_xlim([0,8*atm.ts])
            ax2.set_xlim([0,30000])
            ax2.set_ylabel('Spectral Flux')
            ax2.set_xlabel('Wavenumber')
            ax2.set_title('Spectral OLR')
            # plt.show()
            plt.savefig(output_dir+'/TP_profile_'+str(round(time_current))+'.png', bbox_inches="tight")
            plt.show()
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

    return atm.LW_flux_up[-1]
  
# Dry adjustment routine
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

# Moist adjustment routine
def moistAdj(atm):

    T = atm.temp
    p = atm.p

    if Moist_Adjustment: # switch
        
        Tdew=ga.General_moist_adiabat

        # Downward pass
        for i in range(len(T)-1):
            if T[i]<Tdew(p[i]):
                T[i]=Tdew(p[i])   # Temperature stays the same during phase change
        # Upward pass
        for i in range(len(T)-2,-1,-1): 
            if T[i]<Tdew(p[i]):
                T[i]=Tdew(p[i])

# Time integration for n steps
def steps(atm, stellar_toa_heating):
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
    
    # Dry adjustment step
    for iadj in range(10):
        dryAdj(atm)
        moistAdj(atm)
    Tad = atm.temp[-1]*(atm.p/atm.p[-1])**atm.Rcp
    
    # ** Temporary kludge to keep stratosphere from getting too cold
    atm.temp = np.where(atm.temp<50.,50.,atm.temp)  #**KLUDGE

    # Dummies for separate LW and stellar. **FIX THIS**
    fluxStellar = fluxLW = heatStellar = heatLW = np.zeros(atm.nlev)
    
    return atm
