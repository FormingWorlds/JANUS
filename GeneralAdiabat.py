'''
Created 01/12/19

@authors: 
Ryan Boukrouche (RB)
Tim Lichtenberg (TL)

This file builds a self-consistent temperature profile from Li, Ingersoll & Oyafuso 2018, from the ground up,
using the Runge-Kutta 4 scheme from ClimateUtilities.py.
'''

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
    "CO2_1"   : "#ce5e5e",
    "CO2_2"   : "#B91919",
    "CO2_3"   : "#811111",
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
    "He_3"    : "#465559",
    "NH3_1"   : "#acbbbf",
    "NH3_2"   : "#768E95",
    "NH3_3"   : "#465559"
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

## Dew point temperature [K] given a pressure p [Pa]. Select the molecule of interest with the switch argument (a string).
def Tdew(switch, p): 
    
    # Avoid math error for p = 0
    p = np.max([p, 1e-100])

    if switch == 'H2O':
        Tref = 373.15 # K, boiling point of H2O at 1 atm 
        pref = 1e5 # esat('H2O',Tref) returns 121806.3 Pa, should return 1e5       
        L_H2O=phys.water.L_vaporization*phys.water.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_H2O)*math.log(p/pref))
        #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
    if switch == 'CH4':
        Tref = 148.15 # K, arbitrary point (148.15K,esat(148.15K)=9.66bar) on the L/G coexistence curve of methane 
        pref = p_sat('CH4',Tref)
        L_CH4=phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_CH4)*math.log(p/pref))
    if switch == 'CO2':
        Tref = 253. # K, arbitrary point (253K,esat(253K)=20.9bar) on the coexistence curve of CO2 
        pref = p_sat('CO2',Tref)
        L_CO2=phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_CO2)*math.log(p/pref))
    if switch == 'CO':
        Tref = 100. # K, arbitrary point (100K,esat(100K)=4.6bar) on the coexistence curve of CO 
        pref = p_sat('CO',Tref)
        L_CO=phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_CO)*math.log(p/pref))
    if switch == 'N2':
        Tref = 98.15 # K, arbitrary point (98.15K,esat(98.15K)=7.9bar) on the coexistence curve of N2 
        pref = p_sat('N2',Tref)
        L_N2=phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_N2)*math.log(p/pref))
    if switch == 'O2':
        Tref = 123.15 # K, arbitrary point (123.15K,esat(123.15K)=21.9bar) on the coexistence curve of O2 
        pref = p_sat('O2',Tref)
        L_O2=phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_O2)*math.log(p/pref))
    if switch == 'H2':
        Tref = 23.15 # K, arbitrary point (23.15K,esat(23.15K)=1.7bar) on the coexistence curve of H2 
        pref = p_sat('H2',Tref)
        L_H2=phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_H2)*math.log(p/pref))
    if switch == 'He':
        Tref = 4.22 # K, boiling point of He at 1 atm 
        pref = 1e5 # esat('He',Tref) returns 45196 Pa, should return 1e5
        L_He=phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_He)*math.log(p/pref))
    if switch == 'NH3':
        Tref = 273.15 # K, arbitrary point (273.15K,esat(273.15K)=8.6bar) on the coexistence curve of NH3 
        pref = p_sat('NH3',Tref)
        L_NH3=phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_NH3)*math.log(p/pref))
    
## Molar latent heat [J mol-1] for gas phase considered given a temperature T [K]. 
## Select the molecule of interest with the switch argument (a string).
def L_heat(switch, T, P):

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
        # Numerical issus with 0.?
        L_heat = 0.         #1e-50
        # L_heat = L_vaporization*MolecularWeight*1e-3 

    # No latent heat contribution in moist adiabat if below p_sat
    if P < p_sat(switch, T):
    # if T > Tdew(switch, psat):
        L_heat = 0.

    return L_heat  
    
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

def slopeRay( logpa, logT ):
    eps     = phys.H2O.MolecularWeight/phys.air.MolecularWeight
    L       = phys.H2O.L_vaporization
    Ra      = phys.air.R
    Rc      = phys.H2O.R
    cpa     = phys.air.cp
    cpc     = phys.H2O.cp
    pa      = math.exp(logpa)
    T       = math.exp(logT)
    qsat    = eps*(satvp(T)/pa)
    num     = (1. + (L/(Ra*T))*qsat)*Ra
    den     = cpa + (cpc + (L/(Rc*T) - 1.)*(L/T))*qsat
    return num/den

# # Define pressure levels for SOCRATES (?)
# def set_pressure_array( atm ):
   
#     atm.p     = np.ones(atm.nlev)
#     atm.pl    = np.ones(atm.nlev+1)
#     rat       = (atm.ptop/atm.ps)**(1./atm.nlev)
#     logLevels = [atm.ps*rat**i for i in range(atm.nlev+1)]
#     levels    = [atm.ptop + i*(atm.ps-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
#     atm.pl    = np.array(logLevels)
#     atm.p     = (atm.pl[1:] + atm.pl[:-1]) / 2
    
#     return atm

# Dry adiabat as a comparison
def dry_adiabat( T_surf, p_array, cp_array ):

    T_dry   = np.zeros(len(p_array))
    P_surf  = np.amax(p_array)

    # Check if scalar or cp array
    if not isinstance(cp_array, list):
        cp_array = np.ones(len(p_array))*cp_array

    for idx, prs in enumerate(p_array):

        if cp_array[idx] > 0.:

            T_dry[idx] = T_surf * ( prs / P_surf ) ** ( R_universal / cp_array[idx] )

    return T_dry


# dlnT/dlnP slope function from Li, Ingersoll & 
def moist_slope(lnP, lnT, atm):
    
    # T instead lnT
    tmp = math.exp(lnT)

    # # Kludge for H2 p_sat > 7200 for T > 13.95K
    # if tmp < 13.95: 
    #     tmp = 13.95

    # Find current atm index
    idx = int(np.amax(atm.ifatm))

    # Sum terms in equation
    num_sum     = 0.
    denom_sum1  = 0. 
    denom_sum2  = 0. 
    denom_sum3  = 0.

    # Calculate sums over volatiles
    for vol in atm.vol_list.keys(): 

        # print(vol, atm.x_moist[vol][idx], atm.xd[idx])

        # Coefficients
        eta_vol     = atm.x_gas[vol][idx] / atm.xd[idx]
        beta_vol    = L_heat(vol, tmp, atm.p_vol[vol][idx]) / (R_universal * tmp) 

        # Sum in numerator
        num_sum     += eta_vol * beta_vol

        # Sums in denominator
        denom_sum1  += eta_vol * (beta_vol**2.)
        denom_sum3  += eta_vol
                              
    # Sum 2 in denominator  
    denom_sum2  = num_sum ** 2.

    # Collect terms
    numerator   = 1. + num_sum
    denominator = (atm.cp[idx] / R_universal) + (denom_sum1 + denom_sum2) / (1. + denom_sum3)

    # dlnT/dlnP
    dlnTdlnP = numerator / denominator

    # Moist adiabat slope
    return dlnTdlnP

# Apply condensation and renormalize volatile abundances in gas and condensed phases
def condensation( atm ):

    # Find current level
    idx = int(np.amax(atm.ifatm))

    # Temperature floor
    tmp = np.amax([atm.tmp[idx], 20])

    # Recalculated total pressure
    P_tot_new = 0.

    # Update mixing ratios and partial pressures
    for vol in atm.vol_list.keys():

        # Scaled partial pressure
        atm.p_vol[vol][idx] = atm.vol_list[vol] * atm.p[idx]

        # Saturation vapor pressure
        p_vol_sat           = p_sat(vol, tmp)

        # Condensation if p_old > p_sat: moist species
        if atm.p_vol[vol][idx] > p_vol_sat:

            # Condensate phase molar concentration
            atm.x_cond[vol][idx] = (atm.p_vol[vol][idx] - p_vol_sat) / atm.p[idx]

            # Reduce gas phase molar concentration due to condensation
            atm.x_gas[vol][idx]  = p_vol_sat / atm.p[idx]
            
            # Set species partial pressure to p_sat
            atm.p_vol[vol][idx]  = p_vol_sat

            # Add to molar concentration of condensed phase
            atm.xc[idx]          += atm.x_cond[vol][idx]

            # Add to molar concentration of gas phase
            atm.xv[idx]          += atm.x_gas[vol][idx]

        # Does not condense: dry species
        else:

            # No condensates
            atm.x_cond[vol][idx] = 0.

            # Gas phase molar concentration unchanged
            atm.x_gas[vol][idx]  = atm.p_vol[vol][idx] / atm.p[idx]

            # Add to molar concentration of dry species
            atm.xd[idx]          += atm.x_gas[vol][idx]

        # Update cp w/ molar concentration
        atm.cp[idx]         += (atm.x_gas[vol][idx] + atm.x_cond[vol][idx]) * cpv(vol)
        # atm.cp[idx]         += atm.x_gas[vol][idx] * cpv(vol)

        # Update total pressure
        P_tot_new           += atm.p_vol[vol][idx]

    # # Reset total pressure due to condensation effects
    # atm.p[idx] = P_tot_new

    # Reset surface total pressure at surface
    if idx == 0:
        atm.p[idx] = P_tot_new
        atm.ps     = P_tot_new
        atm.ptop   = np.amin([atm.ps*1e-10, 1e-5])

    # Renormalize cp w/ molar concentration
    atm.cp[idx]  = atm.cp[idx] / (atm.xd[idx] + atm.xv[idx] + atm.xc[idx])

    # Dry concentration floor
    atm.xd[idx] = np.amax([atm.xd[idx], 1e-10])

    ## 'Molar abundance in one mole of heterogeneous gas mixture' (Li, Ingersoll, Oyafuso 2018)

    # Phase abundances
    atm.mrd[idx] = atm.xd[idx] / ( atm.xd[idx] + atm.xv[idx] )
    atm.mrv[idx] = atm.xv[idx] / ( atm.xd[idx] + atm.xv[idx] )
    atm.mrc[idx] = atm.xc[idx] / ( atm.xd[idx] + atm.xv[idx] )

    # Individual volatile abundances
    for vol in atm.vol_list.keys():
        atm.mr_gas[vol][idx]  = atm.x_gas[vol][idx] / ( atm.xd[idx] + atm.xv[idx] )
        atm.mr_cond[vol][idx] = atm.x_cond[vol][idx] / ( atm.xd[idx] + atm.xv[idx] )
        
        # Update cp w/ molar abundance
        atm.cp_mr[idx]         += (atm.mr_gas[vol][idx] + atm.mr_cond[vol][idx]) * cpv(vol)

        # print(idx, vol, atm.x_gas[vol][idx]/atm.xd[idx])

    # Renormalize cp w/ molar abundance
    atm.cp_mr[idx]  = atm.cp[idx] / (atm.mrd[idx] + atm.mrv[idx])

    return atm

# Builds the generalized moist adiabat from slope(lnP, lnT, atm object)
def general_adiabat( atm ):

    ### Initialization

    # Initialize the tuple solution
    moist_tuple    = [] #[tuple([np.log(atm.ps), atm.ts])] 
   
    # Negative increment to go from ps to ptop < ps       
    step            = -.1

    # Integration counter
    idx             = 0  

    # Calculate condensation
    atm             = condensation(atm)

    # Create the integrator instance                                              
    int_slope       = integrator(moist_slope, np.log(atm.ps), np.log(atm.ts), step)

    # Update parameters used in the slope function dlntT/dlnP
    int_slope.setParams(atm)

    # # Integrator for comparison with Ray's slope
    # int_slopeRay    = integrator(slopeRay, np.log(atm.ps), np.log(atm.ts), step)

    ### Integration

    while atm.p[idx] > atm.ptop:
        
        # Execute the Runge-Kutta integrator, fill array of tuples
        moist_tuple.append(int_slope.next())

        # Fill new T,P values
        atm.p[idx+1]    = np.exp(int_slope.x)
        atm.tmp[idx+1]  = np.exp(int_slope.y)

        # print("RK4 step, idx:", idx, round(atm.p[idx+1],5), round(atm.tmp[idx+1],5))

        # Set next level to calculate
        idx             += 1
        atm.ifatm[idx]  = idx

        # Calculate condensation at next level
        atm             = condensation(atm)

        # # Update parameters used in the slope function dT/dlnP
        # int_slope.setParams(atm)

    # Interpolate staggered nodes
    atm.pl      = (atm.p[1:] + atm.p[:-1]) / 2.
    atm.tmpl    = np.interp(atm.pl, np.flip(atm.p), np.flip(atm.tmp))
    for vol in atm.vol_list.keys():
        # atm.x_gasl[vol] = np.zeros(len(atm.pl))
        atm.x_gasl[vol] = np.interp(atm.pl, np.flip(atm.p), np.flip(atm.x_gas[vol]))  

    return atm 


# Plotting
def plot_adiabats(atm):

    sns.set_style("ticks")
    sns.despine()

    ls_moist    = 2.5
    ls_dry      = 2.0
    ls_ind      = 1.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

    # # Ray's moist adiabat as comparison
    # # (1 condensible + 1 non-condensible gas, fixed partial pressure of the non-condensible)
    # ray_vol_cond        = "H2O"
    # ray_vol_noncond     = "N2"
    # moist_adiabat_ray   = phys.MoistAdiabat(phys.H2O,phys.N2)
    # p_dry_gas_surf      = atm.p_vol[ray_vol_noncond][0]
    # print("Ray adiabat settings (P_surf_noncond, T_surf, P_top):", p_dry_gas_surf, atm.ts, atm.ptop)
    # p_ray, T_ray, molarCon_ray, massCon_ray = moist_adiabat_ray( p_dry_gas_surf, atm.ts, atm.ptop )
    # p_ray_interp, T_ray_interp, molarCon_ray_interp, massCon_ray_interp = moist_adiabat_ray(p_dry_gas_surf, atm.ts, atm.ptop, atm.p)
    # # Plot Ray's H2O moist adiabat function
    # ax1.semilogy(T_ray, p_ray, lw=ls_dry, color=vol_colors[ray_vol_cond+"_3"], ls="-.", label=vol_latex[ray_vol_cond]+"/"+vol_latex[ray_vol_cond]+r' adiabat Ray') # label=r'p$_{non-cond.}$ = '+"{:.2f}".format(p_dry_gas)+' Pa'
    # ax2.semilogy(molarCon_ray, p_ray, lw=ls_dry, color=vol_colors[ray_vol_cond+"_3"], ls="-.", label=r"Ray's function")
    
    # For reference p_sat lines
    T_sat_array    = np.linspace(20,3000,1000) 
    p_partial_sum  = np.zeros(len(atm.tmp))

    vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}

    # Individual species
    for vol in vol_list_sorted.keys():

        # Only if volatile is present
        if atm.vol_list[vol] > 0.:
    
            # Saturation vapor pressure for given temperature
            Psat_array = [ p_sat(vol, T) for T in T_sat_array ]
            ax1.semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$'+vol_latex[vol], lw=ls_ind, ls=":", color=vol_colors[vol+"_1"])

            # Plot partial pressures
            ax1.semilogy(atm.tmp, atm.p_vol[vol], color=vol_colors[vol+"_1"], lw=ls_ind, ls="-", label=r'$p$'+vol_latex[vol],alpha=0.99)

            # Sum up partial pressures
            p_partial_sum += atm.p_vol[vol]

            # Plot individual molar concentrations
            ax2.semilogy(atm.x_cond[vol],atm.p, color=vol_colors[vol+"_2"], lw=ls_ind, ls="--", label=vol_latex[vol]+" cloud")
            ax2.semilogy(atm.x_gas[vol],atm.p, color=vol_colors[vol+"_1"], lw=ls_ind, ls="-", label=vol_latex[vol]+" gas")
            
    # Plot sum of partial pressures as check
    ax1.semilogy(atm.tmp, p_partial_sum, color="green", lw=ls_dry, ls="-", label=r'$\sum p_\mathrm{vol}$',alpha=0.99)

    # Dry adiabat as comparison
    ax1.semilogy( dry_adiabat( atm.ts, atm.p, atm.cp ), atm.p , color=vol_colors["black_2"], ls="--", lw=ls_dry, alpha=0.5)              # Condensed abundances
    ax1.semilogy( dry_adiabat( atm.ts, atm.p, atm.cp[0] ), atm.p , color=vol_colors["black_2"], ls="--", lw=ls_dry, label=r'Dry adiabat')   # Initial abundances

    # General moist adiabat
    ax1.semilogy(atm.tmp, atm.p, color=vol_colors["black_1"], lw=ls_moist,label="Moist adiabat",alpha=0.99)

    # Phase molar concentrations
    ax2.semilogy(atm.xd+atm.xv,atm.p, color=vol_colors["black_2"], lw=ls_ind, ls=":", label=r"$X_\mathrm{gas}$")

    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ (K)')
    ax1.set_ylabel(r'Total pressure $P$ (Pa)')
    ax1.set_title('Adiabats & individual Clausius-Clapeyron slopes')
    ax1.legend(ncol=np.min([len(atm.vol_list)+1,2]))
    ax1.set_xlim([0,np.max(atm.ts)])

    ax2.invert_yaxis()
    ax2.set_title('Phase & individual volatile abundances')
    ax2.set_xlabel(r'Molar concentration $X^{\mathrm{vol}}_{\mathrm{phase}}$')
    ax2.set_ylabel(r'Total pressure $P$ (Pa)')
    # ax2.set_title('Relative abundances with condensation')
    ax2.legend(ncol=2)
    ax2.set_xlim(left=-0.05, right=1.05)

    ax1.set_ylim(top=atm.ptop)
    ax1.set_ylim(bottom=atm.ps)
    ax2.set_ylim(top=atm.ptop)
    ax2.set_ylim(bottom=atm.ps)

    ax2.set_xlim([1e-3, 1.05])

    ax2.set_xscale("log")
    
    # plt.show()

    plt.savefig('./output/general_adiabat.pdf', bbox_inches = 'tight')
    plt.close(fig)  

    return


####################################
##### Stand-alone initial conditions
####################################

# Surface pressure & temperature
P_surf                  = 1e+8         # Pa
T_surf                  = 800.         # K

# Volatile molar concentrations: ! must sum to one !
vol_list = { 
              "H2O" : .3, 
              "CO2" : .3,
              "H2"  : .2, 
              "N2"  : .0,  
              "CH4" : .2, 
              "O2"  : .0, 
              "CO"  : .0, 
              "He"  : .0,
              "NH3" : .0, 
            }

# Create atmosphere object
atm                     = atmos(T_surf, P_surf, vol_list)

# Calculate moist adiabat + condensation
atm                     = general_adiabat(atm)

# Plot adiabat
plot_adiabats(atm)

