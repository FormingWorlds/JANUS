# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 13:17:05 2019

@author: Ryan

This file builds a self-consistent temperature profile from Li, Ingersoll & Oyafuso 2018, from the ground up,
based the Euler scheme T[i] = T[i-1] + dlnT[i-1]. The temperature profile is not initialized.

"""

import time
import numpy as np
import math,phys
from ClimateUtilities import *
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from atmosphere_column import atmos
import seaborn as sns

# Color definitions
vol_colors = {
    "black_1" : "#000000",
    "black_2" : "#323232",
    "black_3" : "#7f7f7f",
    "H2O_1"   : "#8db4cb",
    "H2O_2"   : "#4283A9",
    "H2O_3"   : "#274e65",
    "CO2_1"   : "#811111",
    "CO2_2"   : "#B91919",
    "CO2_3"   : "#ce5e5e",
    "H2_1"    : "#a0d2cb",
    "H2_2"    : "#62B4A9",
    "H2_3"    : "#3a6c65",
    "CH4_1"   : "#eb9194",
    "CH4_2"   : "#E6767A",
    "CH4_3"   : "#b85e61",
    "CO_1"    : "#eab597",
    "CO_2"    : "#DD8452",
    "CO_3"    : "#844f31",
    "N2_1"    : "#c29fb2",
    "N2_2"    : "#9A607F",
    "N2_3"    : "#4d303f",  
    "S_1"     : "#f1ca70",
    "S_2"     : "#EBB434",
    "S_3"     : "#a47d24",    
    "O2_1"    : "#57ccda",
    "O2_2"    : "#2EC0D1",
    "O2_3"    : "#2499a7",
    "He_1"    : "#acbbbf",
    "He_2"    : "#768E95",
    "He_3"    : "#465559"
}

# Volatile Latex names
vol_latex = {
    "H2O"   : r"H$_2$O",
    "CO2"   : r"CO$_2$",
    "H2"    : r"H$_2$" ,
    "CH4"   : r"CH$_4$",
    "CO"    : r"CO" ,
    "N2"    : r"N$_2$" ,
    "S"     : r"S"  ,
    "O2"    : r"O$_2$" ,
    "He"    : r"He" ,
    "NH3"   : r"NH$_3$"
}

#--------- Importing thermodynamical properties of gases -----------

R_universal = 8.31446261815324 # Universal gas constant, J.K-1.mol-1, should probably go elsewhere since it's a constant    
satvp = phys.satvps_function(phys.H2O)

## Saturation vapor pressure [Pa] for given temperature T [K]. 
## Assuming the ideal gas law and a constant latent heat of vaporization. 
## Select the molecule of interest with the switch argument (a string).
def p_sat(switch,T): 
    
    # Define volatile
    if switch == 'H2O':
        e = phys.satvps_function(phys.water)
    if switch == 'CH4':
        e = phys.satvps_function(phys.methane)
    if switch == 'CO2':
        e = phys.satvps_function(phys.co2)
    if switch == 'CO':
        e = phys.satvps_function(phys.co)
    if switch == 'N2':
        e = phys.satvps_function(phys.n2)
    if switch == 'O2':
        e = phys.satvps_function(phys.o2)
    if switch == 'H2':
        e = phys.satvps_function(phys.h2)
    if switch == 'He':
        e = phys.satvps_function(phys.he)
    if switch == 'NH3':
        e = phys.satvps_function(phys.nh3)
    
    # Return saturation vapor pressure
    return e(T)
    
#-----------------------------Antoine coefficients-----------------------------
# Need to define what happens when T<Lower_Bound. At p[23]=145Pa we have Tgrid[23,0]=254.6K<Lower_Bound so Tdew[i<23] is not defined
"""
Upper_Bound=573.   # validity range in K from https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
Lower_Bound=255.9   
 
def A(switch,p):
    if switch == 'H2O':        
        if p>=p_sat(Lower_Bound) and p<=p_sat(373.):
            return 4.6543
        if p>=p_sat(379.) and p<=p_sat(Upper_Bound):
            return 3.55959
        
        #if pp>p_sat(373.) and pp<p_sat(379.): # missing values between T=373K and 379K
            #return interp()
        

def B(switch,p):
    if pp>=p_sat(Lower_Bound) and pp<=p_sat(373.):
        return 1435.264
    if pp>=p_sat(379.) and pp<=p_sat(Upper_Bound):
        return 643.748
    
    #if pp>p_sat(373.) and pp<p_sat(379.): # missing values between T=373K and 379K
        #return interp()
    

def C(switch,p):
    if pp>=p_sat(Lower_Bound) and pp<=p_sat(373.):
        return -64.848
    if pp>=p_sat(379.) and pp<=p_sat(Upper_Bound):
        return -198.043
    
    #if pp>p_sat(373.) and pp<p_sat(379.): # missing values between T=373K and 379K
        #return interp()
    
    
    #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
"""  

## Molar latent heat [J mol-1] for gas phase considered given a temperature T [K]. 
## Select the molecule of interest with the switch argument (a string).
def L_heat( switch, T ):

    if switch == 'H2O':
        L_sublimation   = phys.H2O.L_sublimation
        L_vaporization  = phys.H2O.L_vaporization
        MolecularWeight = phys.H2O.MolecularWeight
        T_triple        = phys.H2O.TriplePointT
        T_crit          = phys.H2O.CriticalPointT
    
    if switch == 'CH4':
        L_sublimation   = phys.CH4.L_sublimation
        L_vaporization  = phys.CH4.L_vaporization
        MolecularWeight = phys.CH4.MolecularWeight
        T_triple        = phys.CH4.TriplePointT
        T_crit          = phys.CH4.CriticalPointT
        
    if switch == 'CO2':
        L_sublimation   = phys.CO2.L_sublimation
        L_vaporization  = phys.CO2.L_vaporization
        MolecularWeight = phys.CO2.MolecularWeight
        T_triple        = phys.CO2.TriplePointT
        T_crit          = phys.CO2.CriticalPointT
        
    if switch == 'CO':
        L_sublimation   = phys.CO.L_sublimation
        L_vaporization  = phys.CO.L_vaporization
        MolecularWeight = phys.CO.MolecularWeight
        T_triple        = phys.CO.TriplePointT
        T_crit          = phys.CO.CriticalPointT
            
    if switch == 'N2':
        L_sublimation   = phys.N2.L_sublimation
        L_vaporization  = phys.N2.L_vaporization
        MolecularWeight = phys.N2.MolecularWeight
        T_triple        = phys.N2.TriplePointT
        T_crit          = phys.N2.CriticalPointT
        
    if switch == 'O2':
        L_sublimation   = phys.O2.L_sublimation
        L_vaporization  = phys.O2.L_vaporization
        MolecularWeight = phys.O2.MolecularWeight
        T_triple        = phys.O2.TriplePointT
        T_crit          = phys.O2.CriticalPointT
    
    if switch == 'H2':
        L_sublimation   = phys.H2.L_vaporization # No H2 sublimation
        L_vaporization  = phys.H2.L_vaporization
        MolecularWeight = phys.H2.MolecularWeight
        T_triple        = phys.H2.TriplePointT
        T_crit          = phys.H2.CriticalPointT
        
    if switch == 'He':
        L_sublimation   = phys.He.L_vaporization  # No He sublimation
        L_vaporization  = phys.He.L_vaporization
        MolecularWeight = phys.He.MolecularWeight
        T_triple        = phys.He.TriplePointT
        T_crit          = phys.He.CriticalPointT
        
    if switch == 'NH3':
        L_sublimation   = phys.NH3.L_sublimation
        L_vaporization  = phys.NH3.L_vaporization
        MolecularWeight = phys.NH3.MolecularWeight
        T_triple        = phys.NH3.TriplePointT
        T_crit          = phys.NH3.CriticalPointT



    # Gas-solid transition
    if T <= T_triple:
        # Conversion from J.kg-1 to J.mol-1 (molecular weight is in g/mol in phys.f90)
        L_heat = L_sublimation*MolecularWeight*1e-3 
    # Gas-liquid transition
    elif T <= T_crit:
        L_heat = L_vaporization*MolecularWeight*1e-3 
    # Super-critical state
    else:
        # Numerical issus with 0.
        L_heat = 1e-30

    return L_heat  
    
## Dew point temperature [K] array for given pressure array [Pa] and surface T [K]. 
## Select the molecule of interest with the switch argument (a string)
def Tdew(switch, prs, T_surf, L_use): 

    # Calculate dew-point for each pressure
    T_dew = T_surf / ( 1. - ( T_surf * R_universal / L_use ) * np.log( prs / p_sat(switch,T_surf) ) )

    # # Re-calc with L_heat switch
    # for idx, T in enumerate(T_dew):
    #     L_use = L_heat( switch, T )
    #     T_dew[idx] = T_surf / ( 1. - ( T_surf * R_universal / L_use ) * np.log( prs[idx] / p_sat(switch,T_surf) ) )

    return T_dew
    
def cpv(switch): 
    """Molar heat capacities [J.K-1.mol-1]. Select the molecule of interest with the switch argument (a string)."""
    if switch == 'H2O':
        return phys.water.cp*phys.water.MolecularWeight*1e-3
    if switch == 'CH4':
        return phys.methane.cp*phys.methane.MolecularWeight*1e-3
    if switch == 'CO2':
        return phys.co2.cp*phys.co2.MolecularWeight*1e-3
    if switch == 'CO':
        return phys.co.cp*phys.co.MolecularWeight*1e-3
    if switch == 'N2':
        return phys.n2.cp*phys.n2.MolecularWeight*1e-3
    if switch == 'O2':
        return phys.o2.cp*phys.o2.MolecularWeight*1e-3
    if switch == 'H2':
        return phys.h2.cp*phys.h2.MolecularWeight*1e-3
    if switch == 'He':
        return phys.he.cp*phys.he.MolecularWeight*1e-3
    if switch == 'NH3':
        return phys.nh3.cp*phys.nh3.MolecularWeight*1e-3    
    
def slope(lnP,lnT,params):
    """Returns the slope dT/dlnP given a temperature T [K], the natural logarithm of a pressure [Pa] lnP and the dictionary atm_chemistry_arrays passed to this function with the ClimateUtilities object params."""
    
    T = math.exp(lnT)
    # Kludge for H2 esat>7200 for T>13.95K
    if T < 13.95: 
        T = 13.95
    
    # Define individual abundances from the last abundance calculated
    xH2O = params.atm_chemistry_arrays['H2O'][-1]
    xCO2 = params.atm_chemistry_arrays['CO2'][-1]
    xCH4 = params.atm_chemistry_arrays['CH4'][-1]
    xCO  = params.atm_chemistry_arrays['CO'][-1]
    xN2  = params.atm_chemistry_arrays['N2'][-1]
    xO2  = params.atm_chemistry_arrays['O2'][-1]
    xH2  = params.atm_chemistry_arrays['H2'][-1]
    xHe  = params.atm_chemistry_arrays['He'][-1]
    xNH3 = params.atm_chemistry_arrays['NH3'][-1]
    xd   = 1. - ( xH2O + xCO2 + xCH4 + xCO + xN2 + xO2 + xH2 + xHe + xNH3 )
    
    # Molar heat capacity of the dry species
    cpd = phys.air.cp*phys.air.MolecularWeight*1.e-3 #cpv('N2')

    # Avoid division by zero if atmosphere consists of condensibles only
    xd = np.max([xd,1e-8])
                         
    xv_cpv = xH2O * cpv('H2O') + \
             xCO2 * cpv('CO2') + \
             xCH4 * cpv('CH4') + \
             xCO  * cpv('CO' ) + \
             xN2  * cpv('N2' ) + \
             xO2  * cpv('O2' ) + \
             xH2  * cpv('H2' ) + \
             xHe  * cpv('He' ) + \
             xNH3 * cpv('NH3')
  
    first_term = (xH2O/xd) * ( L_heat('H2O',T) / (R_universal*T) )**2 + \
                 (xCO2/xd) * ( L_heat('CO2',T) / (R_universal*T) )**2 + \
                 (xCH4/xd) * ( L_heat('CH4',T) / (R_universal*T) )**2 + \
                 (xCO /xd) * ( L_heat('CO' ,T) / (R_universal*T) )**2 + \
                 (xN2 /xd) * ( L_heat('N2' ,T) / (R_universal*T) )**2 + \
                 (xO2 /xd) * ( L_heat('O2' ,T) / (R_universal*T) )**2 + \
                 (xH2 /xd) * ( L_heat('H2' ,T) / (R_universal*T) )**2 + \
                 (xHe /xd) * ( L_heat('He' ,T) / (R_universal*T) )**2 + \
                 (xNH3/xd) * ( L_heat('NH3',T) / (R_universal*T) )**2
    
    quadr_term=(
                 (xH2O/xd) * L_heat('H2O',T) / (R_universal*T) + \
                 (xCO2/xd) * L_heat('CO2',T) / (R_universal*T) + \
                 (xCH4/xd) * L_heat('CH4',T) / (R_universal*T) + \
                 (xCO /xd) * L_heat('CO' ,T) / (R_universal*T) + \
                 (xN2 /xd) * L_heat('N2' ,T) / (R_universal*T) + \
                 (xO2 /xd) * L_heat('O2' ,T) / (R_universal*T) + \
                 (xH2 /xd) * L_heat('H2' ,T) / (R_universal*T) + \
                 (xHe /xd) * L_heat('He' ,T) / (R_universal*T) + \
                 (xNH3/xd) * L_heat('NH3',T) / (R_universal*T)             
                )**2
        
    sum_abundances = xH2O + xCO2 + xCH4 + xCO + xN2 + xO2 + xH2 + xHe + xNH3
                                             
    sum_ratio_abundances = xH2O / xd + xCO2 / xd + \
                                       xCH4 / xd + \
                                       xCO  / xd + \
                                       xN2  / xd + \
                                       xO2  / xd + \
                                       xH2  / xd + \
                                       xHe  / xd + \
                                       xNH3 / xd
                                                      
    num_sum = (xH2O/xd) * L_heat('H2O',T) / (R_universal*T) + \
              (xCO2/xd) * L_heat('CO2',T) / (R_universal*T) + \
              (xCH4/xd) * L_heat('CH4',T) / (R_universal*T) + \
              (xCO /xd) * L_heat('CO' ,T) / (R_universal*T) + \
              (xN2 /xd) * L_heat('N2' ,T) / (R_universal*T) + \
              (xO2 /xd) * L_heat('O2' ,T) / (R_universal*T) + \
              (xH2 /xd) * L_heat('H2' ,T) / (R_universal*T) + \
              (xHe /xd) * L_heat('He' ,T) / (R_universal*T) + \
              (xNH3/xd) * L_heat('NH3',T) / (R_universal*T)
     
    num       = (1.+ num_sum)
    denom     = (1./R_universal) * (xd*cpd+xv_cpv) / (xd+sum_abundances) + \
                (first_term + quadr_term) / (1.+sum_ratio_abundances)
    
    params.general_cp = (xd*cpd + xv_cpv) / (xd + sum_abundances)
    
    return num/denom  

def slopeRay(logpa,logT):
    eps = phys.H2O.MolecularWeight/phys.air.MolecularWeight
    L = phys.H2O.L_vaporization
    Ra = phys.air.R
    Rc = phys.H2O.R
    cpa = phys.air.cp
    cpc = phys.H2O.cp
    pa = math.exp(logpa)
    T = math.exp(logT)
    qsat = eps*(satvp(T)/pa)
    num = (1. + (L/(Ra*T))*qsat)*Ra
    den = cpa + (cpc + (L/(Rc*T) - 1.)*(L/T))*qsat
    return num/den

# Define pressure levels
def set_pressure_array(atm):
   
    atm.p     = np.ones(atm.nlev)
    atm.pl    = np.ones(atm.nlev+1)
    rat       = (atm.ptop/atm.ps)**(1./atm.nlev)
    logLevels = [atm.ps*rat**i for i in range(atm.nlev+1)]
    levels    = [atm.ptop + i*(atm.ps-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl    = np.array(logLevels)
    atm.p     = (atm.pl[1:] + atm.pl[:-1]) / 2
    
    return atm

def dry_adiabat( T_surf, p_array, Rcp ):

    T_p = T_surf * np.power( ( p_array / np.max(p_array) ), Rcp )

    return T_p

def condensation( atm ):

    # Find new level to calculate
    idx = np.amax(atm.ifatm)

    print("Condensation (lvl ", idx, ")", end=" ")

    # Figure out which species are condensing
    for vol in atm.vol_list:

        print(vol, end=" ")
        
        # Current partial pressure & mixing ratios
        p_vol_old      = atm.p_vol[vol][idx]
        x_moist_old    = atm.x_moist[vol][idx]

        # Saturation vapor pressure at current temperature
        p_vol_sat      = p_sat(vol, atm.tmp[idx])

        # Condensation if p_old > p_sat
        if p_vol_old > p_vol_sat:

            # Reduce species partial pressure to p_sat
            atm.p_vol[vol][idx]   = p_vol_sat

            # Recalculate mixing ratios
            atm.x_moist[vol][idx] = p_vol_sat / P_tot
            atm.x_cond[vol][idx]  = x_moist_old - atm.x_moist[vol][idx]
            atm.x_dry[vol][idx]   = 0.

            # Flag as condensing species
            atm.vol_cond[idx].append(vol)

        # Does not condense: dry species
        else:

            # Add partial pressure to combined partial pressure of dry species
            atm.xd[idx]           += p_vol_old

            # Recalculate mixing ratios
            atm.x_moist[vol][idx] = 0.
            atm.x_cond[vol][idx]  = 0.
            atm.x_dry[vol][idx]   = p_vol_old / atm.p[idx]

            # Flag as dry species
            atm.vol_dry[idx].append(vol)

    # ?? Recalc total pressure ??
    # P_tot = np.sum(atm.p_vol.values())

    print()

    return atm
                        
def moist_adiabat(atm):
    """Builds the generalized moist adiabat from slope(lnP,T,params) and plots the adiabat along with the relative abundances."""   

    # for mode in [ "original" ]: # "hybrid", , "ray1"

    # params = Dummy() # initialize parameter object  

    # print(params)

    # params.atm_chemistry_arrays = {}
    # for x in atm_chemistry:
    #     params.atm_chemistry_arrays['%s' % x] = [] 
          
    # Initialization

    # Initialize the tuple solution
    moist_w_cond    = [] #[tuple([np.log(atm.ps), atm.ts])] 

    # Initialize the final pressure array   
    pL              = []                                                  
    psatL           = []

    # Initialize the exponent of the dry adiabat
    Rcp             = []    

    # Initialize the final temperature array                                             
    TL              = []
   
    # Negative increment to go from ps to ptop < ps       
    step            = -.1

    # Calculate condensation
    atm             = condensation(atm)

    # Create the integrator instance                                              
    int_slope       = integrator(slope, np.log(atm.ps), np.log(atm.ts), step)

    # Integrator for comparison with Ray's slope
    int_slopeRay    = integrator(slopeRay, np.log(atm.ps), np.log(atm.ts), step)

    # Update parameters used in the slope function dT/dlnP
    int_slope.setParams(atm)

    # Integration counter
    idx           = 0    

    logT            = int_slope.y
    logP            = int_slope.x
    Temperature     = np.exp(logT)
    Pressure        = np.exp(logP)
  
    ##### Integration
    while Pressure > atm.ptop:
        
        # Loop on the volatiles 
        for vol in atm_chemistry:                       
            
            # Initialize volatile from original abundance
            params.atm_chemistry_arrays[vol].append(atm_chemistry[vol]) 

            # Original partial pressure
            p_vol = params.atm_chemistry_arrays[vol][-1]*Pressure

            # Saturation vapor pressure for current level
            p_vol_sat  = p_sat(vol,Temperature)

            # Reset partial pressures to saturation vapor pressure
            # H2 is background dry gas, assumed to not condense
            if vol != "H2":
                params.atm_chemistry_arrays[vol][-1] = np.min([p_vol,p_vol_sat])/Pressure

            # if vol != "H2":
                

            # If int_slope.y < Tdew, how can numpy.exp(int_slope.x) < p_sat?
             # min to avoid p_sat > p_molecule       
            # params.atm_chemistry_arrays[vol][-1] = p_sat/Pressure
        
        # Execute the Runge-Kutta integrator, fill array of tuples
        moist_w_cond.append(int_slope.next()) 
        pL.append(numpy.exp(int_slope.x))
        TL.append(numpy.exp(int_slope.y))

        logT        = int_slope.y
        logP        = int_slope.x
        Temperature = np.exp(logT)
        Pressure    = np.exp(logP)
        Rcp.append(R_universal/params.general_cp)

        idx += 1


    # Convert to numpy array. 
    moist_w_cond = numpy.array(moist_w_cond) # lnP accessed through moist_w_cond[:,0]
                                             # T   accessed through moist_w_cond[:,1] 

    pL                 = numpy.array(pL)
    TL                 = numpy.array(TL)    
    abundance_arrayH2O = numpy.array(params.atm_chemistry_arrays['H2O'])
    abundance_arrayNH3 = numpy.array(params.atm_chemistry_arrays['NH3'])
    abundance_arrayCO2 = numpy.array(params.atm_chemistry_arrays['CO2'])
    abundance_arrayCH4 = numpy.array(params.atm_chemistry_arrays['CH4'])
    abundance_arrayCO  = numpy.array(params.atm_chemistry_arrays['CO'])
    abundance_arrayN2  = numpy.array(params.atm_chemistry_arrays['N2'])
    abundance_arrayO2  = numpy.array(params.atm_chemistry_arrays['O2'])
    abundance_arrayH2  = numpy.array(params.atm_chemistry_arrays['H2'])
    abundance_arrayHe  = numpy.array(params.atm_chemistry_arrays['He'])

    

    # Ray's moist adiabat (1 condensible + 1 non-condensible gas, fixed partial pressure of the non-condensible)
    ray_vol_cond = "H2O"
    moist_adiabat_ray = phys.MoistAdiabat(phys.H2O,phys.N2)
    p_dry_gas         = atm.ps*(1.-atm_chemistry[ray_vol_cond])
    p_ray, T_ray, molarCon_ray, massCon_ray = moist_adiabat_ray( p_dry_gas, atm.ts, atm.ptop )
    p_ray_interp, T_ray_interp, molarCon_ray_interp, massCon_ray_interp = moist_adiabat_ray(p_dry_gas, atm.ts, atm.ptop, p_plot)



    ### Plotting

    # Pressure array for interpolation scheme
    p_plot = np.exp(np.linspace(np.log(atm.ptop),np.log(atm.ps),100))
    sns.set_style("ticks")
    sns.despine()

    ls_adiabat  = 2
    ls_ind      = 1.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

    # Plot Ray's H2O moist adiabat function
    ax1.semilogy(T_ray, p_ray, lw=ls_adiabat, color=vol_colors[ray_vol_cond+"_2"], ls="--", label=vol_latex[ray_vol_cond]+r' adiabat Ray') # label=r'p$_{non-cond.}$ = '+"{:.2f}".format(p_dry_gas)+' Pa'
    ax2.semilogy(molarCon_ray, p_ray, lw=ls_adiabat, color=vol_colors[ray_vol_cond+"_2"], ls="--", label=r"Ray's function")
    
    
    # Interpolate the integrated adiabat to p_plot?
    Interpolate = False
    
    if Interpolate:        
        #pL=numpy.flip(pL)
        T1                  = interp(pL,TL)
        abundance_arrayH2O1 = interp(pL,abundance_arrayH2O)
        abundance_arrayNH31 = interp(pL,abundance_arrayNH3)
        abundance_arrayCO21 = interp(pL,abundance_arrayCO2)
        abundance_arrayCH41 = interp(pL,abundance_arrayCH4)
        abundance_arrayCO1  = interp(pL,abundance_arrayCO)
        abundance_arrayN21  = interp(pL,abundance_arrayN2)
        abundance_arrayO21  = interp(pL,abundance_arrayO2)
        abundance_arrayH21  = interp(pL,abundance_arrayH2)
        abundance_arrayHe1  = interp(pL,abundance_arrayHe)
        
        TL                 = numpy.array([T1(pp) for pp in pL])
        abundance_arrayH2O = numpy.array([abundance_arrayH2O1(pp) for pp in pL])
        abundance_arrayNH3 = numpy.array([abundance_arrayNH31(pp) for pp in pL])
        abundance_arrayCO2 = numpy.array([abundance_arrayCO21(pp) for pp in pL])
        abundance_arrayCH4 = numpy.array([abundance_arrayCH41(pp) for pp in pL])
        abundance_arrayCO  = numpy.array([abundance_arrayCO1(pp) for pp in pL])
        abundance_arrayN2  = numpy.array([abundance_arrayN21(pp) for pp in pL])
        abundance_arrayO2  = numpy.array([abundance_arrayO21(pp) for pp in pL])
        abundance_arrayH2  = numpy.array([abundance_arrayH21(pp) for pp in pL])
        abundance_arrayHe  = numpy.array([abundance_arrayHe1(pp) for pp in pL])

        pL = p_plot

    # For reference p_sat lines
    T_sat_array    = np.linspace(0,3000,1000) 

    # Plot single-species dew points
    for vol in [ 'H2O' ]:
        
        # # Dew-point temperature for given pressure
        # L_use = L_heat( vol, atm.ts )
        # Tdew_array = Tdew( vol, pL, atm.ts, L_use )
        # ax1.semilogy( Tdew_array, pL, label=vol_latex[vol]+' dew-point', lw=ls_ind, ls="--", color=vol_colors[vol+"_1"])

        # Saturation vapor pressure for given temperature
        Psat_array = [ p_sat('H2O', T) for T in T_sat_array ]
        ax1.semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$ '+vol_latex[vol], lw=ls_ind, ls=":", color=vol_colors[vol+"_1"])


    #ax1.semilogy(moist_w_cond[:,1],p_plot,color="red",lw=ls_adiabat,label=r'Moist adiabat')

    # if mode ==  "original":
    #     color = "red"
    # if mode ==  "hybrid":
    #     color = "blue"
    # if mode ==  "ray1":
    #     color = "green"
    # label = mode
        
    # Dry adiabat as comparison
    T_dry = dry_adiabat( atm.ts, pL, Rcp )
    ax1.semilogy( T_dry, pL , color='black', ls="--", lw=ls_adiabat, label=r'Dry adiabat')
        
    alpha=0.0
    
    
    xplot_H2  = ax2.semilogy(abundance_arrayH2,pL, color=vol_colors["H2_1"] ,label=r'H$_2$ (dry gas)')
    xplot_H2O = ax2.semilogy(abundance_arrayH2O,pL, color=vol_colors["H2O_1"], label=r'H$_2$O')
    xplot_CO2 = ax2.semilogy(abundance_arrayCO2,pL, alpha=alpha)#,label=r'CO$_2$', alpha=alpha)#
    xplot_CH4 = ax2.semilogy(abundance_arrayCH4,pL, alpha=alpha)#,label=r'CH$_4$', alpha=alpha)
    xplot_CO  = ax2.semilogy(abundance_arrayCO,pL, alpha=alpha)#,label=r'CO', alpha=alpha)
    xplot_N2  = ax2.semilogy(abundance_arrayN2,pL, alpha=alpha)#,label=r'N$_2$', alpha=alpha)
    xplot_O2  = ax2.semilogy(abundance_arrayO2,pL, alpha=alpha)#,label=r'O$_2$', alpha=alpha)

    # Plot moist adiabat
    alpha=0.99

    ax1.semilogy(TL, pL, color=vol_colors["CH4_2"], lw=ls_adiabat,label="Moist adiabat",alpha=alpha) # non-interpolated

    #ax2.semilogy(molarCon_interp,p_plot,lw=ls_adiabat, color="blue", ls="--")
    #ax2.semilogy(molarCon,p_ray,lw=ls_adiabat, color="blue", ls="-",label=r'Ray_molarCon')
    
    # Constant abundances w/o condensation as comparison
    
    # if condensation == True:
    
    # if use_vulcan == 0:
    #     const_abundances_H2O = np.ones(len(p_plot))*atm_chemistry["H2O"]
    #     const_abundances_CO2 = np.ones(len(p_plot))*atm_chemistry["CO2"]
    #     const_abundances_H2  = np.ones(len(p_plot))*atm_chemistry["H2"]
    #     const_abundances_CH4 = np.ones(len(p_plot))*atm_chemistry["CH4"]
    #     const_abundances_CO  = np.ones(len(p_plot))*atm_chemistry["CO"]
    #     const_abundances_N2  = np.ones(len(p_plot))*atm_chemistry["N2"]
    #     const_abundances_O2  = np.ones(len(p_plot))*atm_chemistry["O2"]
    #     const_abundances_He  = np.ones(len(p_plot))*atm_chemistry["He"]
    #     const_abundances_NH3 = np.ones(len(p_plot))*atm_chemistry["NH3"]
    """
    ax2.semilogy(const_abundances_H2O, p_plot, ls="--", color=xplot_H2O[0].get_color() )
    ax2.semilogy(const_abundances_CO2, p_plot, ls="--", color=xplot_CO2[0].get_color() )
    ax2.semilogy(const_abundances_H2,  p_plot, ls="--", color=xplot_H2[0].get_color() )
    ax2.semilogy(const_abundances_CH4, p_plot, ls="--", color=xplot_CH4[0].get_color() )
    ax2.semilogy(const_abundances_CO,  p_plot, ls="--", color=xplot_CO[0].get_color() )
    ax2.semilogy(const_abundances_N2,  p_plot, ls="--", color=xplot_N2[0].get_color() )
    ax2.semilogy(const_abundances_O2,  p_plot, ls="--", color=xplot_O2[0].get_color() )
    """

    # Plot surface pressure
    ax1.axhline(atm.ps, ls=":", lw=ls_ind, color=vol_colors["black_3"])
    ax2.axhline(atm.ps, ls=":", lw=ls_ind, color=vol_colors["black_3"])

    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ (K)')
    ax1.set_ylabel(r'Total pressure $P$ (Pa)')
    ax1.set_title('Adiabats + individual Clausius-Clapeyron slopes')
    ax1.legend(ncol=1)
    ax1.set_xlim([0,np.max(atm.ts)])
    

    # ax2.semilogy(xHe_array,p_plot,label=r'He')
    # ax2.semilogy(xNH3_array,p_plot,label=r'NH$_3$')
    ax2.invert_yaxis()
    ax2.set_title('Mixing ratios')
    ax2.set_xlabel(r'Molar concentration $X_{\mathrm{cond.}}/X_{total}$')
    ax2.set_ylabel(r'Total pressure $P$ (Pa)')
    # ax2.set_title('Relative abundances with condensation')
    ax2.legend(ncol=1)
    ax2.set_xlim(left=-0.05, right=1.05)

    # ax1.set_ylim([1e5,1e-5])
    ax1.set_ylim(top=1e-5)
    ax1.set_ylim(bottom=1e+7)
    ax2.set_ylim(top=1e-5)
    ax2.set_ylim(bottom=1e+7)
    
    # plt.show()
    plt.savefig('./output/general_adiabat.pdf', bbox_inches = 'tight')
    #plt.close(fig)

    return moist_w_cond   


### Stand-alone initial conditions ###

# Create atmosphere object
atm              = atmos()

# Surface temperature of planet
atm.ts           = 600.           # K
atm.tmp[0]       = atm.ts         # K

# Surface pressure
atm.ps           = 1e+5           # Pa
atm.p[0]         = atm.ps         # Pa

# Top pressure
atm.ptop         = atm.ps*1e-10   # Pa

# # Instantiate pressure array
# atm              = set_pressure_array(atm)

# Initialize temperature profile on dry adiabat with 'Earth air'
atm.temp         = atm.ts*(atm.p/atm.p[0])**atm.Rcp

# Initial gas phase molar mixing ratios at surface
atm.x_moist["H2"]  = [ 0.5 ]
atm.x_moist["N2"]  = [ 0. ]
atm.x_moist["H2O"] = [ 0.5 ]
atm.x_moist["CO2"] = [ 0. ]
atm.x_moist["CH4"] = [ 0. ]
atm.x_moist["O2"]  = [ 0. ]
atm.x_moist["CO"]  = [ 0. ]
atm.x_moist["NH3"] = [ 0. ]
atm.x_moist["He"]  = [ 0. ]

# Initialize level-dependent quantities
atm.vol_list        = atm.x_moist.keys() # List of all species for looping
atm.xd              = [ 0. ]           # Combined molar mixing ratio of all 'dry' gas
atm.vol_dry         = [ [] ]           # Dry species per pressure level
atm.vol_cond        = [ [] ]           # Condensing species per pressure level

# Initialize dry and condensed species
for vol in atm.vol_list:
    atm.x_dry[vol]  = [ 0. ]
    atm.x_cond[vol] = [ 0. ]

# Calculate surface partial pressures
for vol in atm.vol_list:
    atm.p_vol[vol]  = [ atm.ps * atm.x_moist[vol][-1] ]
    atm.p[0]        += atm.p_vol[vol]

# Execut moist adiabat function
moist_w_cond     = moist_adiabat(atm)

