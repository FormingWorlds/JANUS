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

#--------- Importing thermodynamical properties of gases -----------

R_universal = 8.31446261815324 # Universal gas constant, J.K-1.mol-1, should probably go elsewhere since it's a constant    

def esat(switch,T): 
    """Saturation vapor pressure [Pa] given a temperature T [K]. Assuming the ideal gas law and a constant latent heat. Select the molecule of interest with the switch argument (a string)."""
    if switch == 'H2O':
        e=phys.satvps_function(phys.water)
        return e(T) 
    if switch == 'CH4':
        e=phys.satvps_function(phys.methane)
        return e(T)
    if switch == 'CO2':
        e=phys.satvps_function(phys.co2)
        return e(T)
    if switch == 'CO':
        e=phys.satvps_function(phys.co)
        return e(T)
    if switch == 'N2':
        e=phys.satvps_function(phys.n2)
        return e(T)
    if switch == 'O2':
        e=phys.satvps_function(phys.o2)
        return e(T)
    if switch == 'H2':
        e=phys.satvps_function(phys.h2)
        return e(T)
    if switch == 'He':
        e=phys.satvps_function(phys.he)
        return e(T)
    if switch == 'NH3':
        e=phys.satvps_function(phys.nh3)
        return e(T)
    
#-----------------------------Antoine coefficients-----------------------------
# Need to define what happens when T<Lower_Bound. At p[23]=145Pa we have Tgrid[23,0]=254.6K<Lower_Bound so Tdew[i<23] is not defined
"""
Upper_Bound=573.   # validity range in K from https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
Lower_Bound=255.9   
 
def A(switch,p):
    if switch == 'H2O':        
        if p>=esat(Lower_Bound) and p<=esat(373.):
            return 4.6543
        if p>=esat(379.) and p<=esat(Upper_Bound):
            return 3.55959
        
        #if pp>esat(373.) and pp<esat(379.): # missing values between T=373K and 379K
            #return interp()
        

def B(switch,p):
    if pp>=esat(Lower_Bound) and pp<=esat(373.):
        return 1435.264
    if pp>=esat(379.) and pp<=esat(Upper_Bound):
        return 643.748
    
    #if pp>esat(373.) and pp<esat(379.): # missing values between T=373K and 379K
        #return interp()
    

def C(switch,p):
    if pp>=esat(Lower_Bound) and pp<=esat(373.):
        return -64.848
    if pp>=esat(379.) and pp<=esat(Upper_Bound):
        return -198.043
    
    #if pp>esat(373.) and pp<esat(379.): # missing values between T=373K and 379K
        #return interp()
    
    
    #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
"""    
    

def Tdew(switch,p): 
    """Dew point temperature [K] given a pressure p [Pa]. Select the molecule of interest with the switch argument (a string)."""
    if switch == 'H2O':
        Tref = 373.15 # K, boiling point of H2O at 1 atm 
        pref = 1e5 # esat('H2O',Tref) returns 121806.3 Pa, should return 1e5       
        L_H2O=phys.water.L_vaporization*phys.water.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_H2O)*math.log(p/pref))
        #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
    if switch == 'CH4':
        Tref = 148.15 # K, arbitrary point (148.15K,esat(148.15K)=9.66bar) on the L/G coexistence curve of methane 
        pref = esat('CH4',Tref)
        L_CH4=phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_CH4)*math.log(p/pref))
    if switch == 'CO2':
        Tref = 253. # K, arbitrary point (253K,esat(253K)=20.9bar) on the coexistence curve of CO2 
        pref = esat('CO2',Tref)
        L_CO2=phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_CO2)*math.log(p/pref))
    if switch == 'CO':
        Tref = 100. # K, arbitrary point (100K,esat(100K)=4.6bar) on the coexistence curve of CO 
        pref = esat('CO',Tref)
        L_CO=phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_CO)*math.log(p/pref))
    if switch == 'N2':
        Tref = 98.15 # K, arbitrary point (98.15K,esat(98.15K)=7.9bar) on the coexistence curve of N2 
        pref = esat('N2',Tref)
        L_N2=phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_N2)*math.log(p/pref))
    if switch == 'O2':
        Tref = 123.15 # K, arbitrary point (123.15K,esat(123.15K)=21.9bar) on the coexistence curve of O2 
        pref = esat('O2',Tref)
        L_O2=phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_O2)*math.log(p/pref))
    if switch == 'H2':
        Tref = 23.15 # K, arbitrary point (23.15K,esat(23.15K)=1.7bar) on the coexistence curve of H2 
        pref = esat('H2',Tref)
        L_H2=phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_H2)*math.log(p/pref))
    if switch == 'He':
        Tref = 4.22 # K, boiling point of He at 1 atm 
        pref = 1e5 # esat('He',Tref) returns 45196 Pa, should return 1e5
        L_He=phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_He)*math.log(p/pref))
    if switch == 'NH3':
        Tref = 273.15 # K, arbitrary point (273.15K,esat(273.15K)=8.6bar) on the coexistence curve of NH3 
        pref = esat('NH3',Tref)
        L_NH3=phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3
        return Tref/(1.-(Tref*R_universal/L_NH3)*math.log(p/pref))

def L_heat(switch,T): 
    """Molar latent heat [J mol-1] for gases considered given a temperature T [K]. Select the molecule of interest with the switch argument (a string)."""
    if switch == 'H2O':
        if T<=Tdew('H2O',esat('H2O',T)):
            return phys.water.L_vaporization*phys.water.MolecularWeight*1e-3 # Conversion from J.kg-1 to J.mol-1 (molecular weight is in g/mol in phys.f90)
        else:
            return 0. # Without condensation, neglect the L terms
    
    if switch == 'CH4':
        if T<=Tdew('CH4',esat('CH4',T)):
            return phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
        else:
            return 0.
        
    if switch == 'CO2':        
        if T<=Tdew('CO2',esat('CO2',T)):
            return phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
        else:
            return 0.
        
    if switch == 'CO':
        if T<=Tdew('CO',esat('CO',T)):
            return phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
        else:
            return 0.
            
    if switch == 'N2':
        if T<=Tdew('N2',esat('N2',T)):
            return phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
        else:
            return 0.
        
    if switch == 'O2':
        if T<=Tdew('O2',esat('O2',T)):
            return phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
        else:
            return 0.
    
    if switch == 'H2':
        if T<=Tdew('H2',esat('H2',T)):
            return phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
        else:
            return 0.
    
    if switch == 'He':
        if T<=Tdew('He',esat('He',T)):
            return phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
        else:
            return 0.
        
    if switch == 'NH3':
        if T<=Tdew('He',esat('He',T)):
            return phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3
        else:
            return 0.
    
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
    
def slope(lnP,T,params):
    """Returns the slope dT/dlnP and builds the abundance arrays given a temperature T [K], the natural logarithm of a pressure [Pa] lnP and the dictionary atm_chemistry_arrays passed to this function with the ClimateUtilities object params."""
    
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
    cpd = cpv('N2') #1000.*phys.air.MolecularWeight*1.e-3 

    # Avoid division by zero if atmosphere consists of condensibles only
    xd = np.max([xd,1e-8])
                         
    xv_cpv = xH2O * cpv('H2O') + xCO2 * cpv('CO2') + \
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
     
    num       = (1.+ num_sum)*T                          # Factor T to get dT/dlnP instead of dlnT/dlnP
    denom     = (1./R_universal) * (xd*cpd+xv_cpv) / (xd+sum_abundances) + \
                (first_term + quadr_term) / (1.+sum_ratio_abundances)
    
    dTdlnP = num/denom    

    return dTdlnP

def set_pressure_array(atm):
    rat       = (atm.ptop/atm.ps)**(1./atm.nlev)
    logLevels = [atm.ps*rat**i for i in range(atm.nlev+1)]
    levels    = [atm.ptop + i*(atm.ps-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl    = np.array(logLevels)
    atm.p     = (atm.pl[1:] + atm.pl[:-1]) / 2
    return atm
                        
def solve_general_adiabat(atm, atm_chemistry, use_vulcan, condensation):
    """Builds the generalized moist adiabat from slope(lnP,T,params) and plots the adiabat along with the relative abundances.""" 
    
    params = Dummy() # initialize parameter object
          
    """
    # Initialize relative surface abundances 
    xH2O = atm_chemistry["H2O"]
    xCO2 = atm_chemistry["CO2"]
    xH2  = atm_chemistry["H2"]
    xCH4 = atm_chemistry["CH4"]
    xCO  = atm_chemistry["CO"]
    xN2  = atm_chemistry["N2"]
    xO2  = atm_chemistry["O2"]
    xHe  = atm_chemistry["He"]
    xNH3 = atm_chemistry["NH3"]
    xd   = 1. - ( xH2O + xCO2 + xCH4 + xCO + xN2 + xO2 + xH2 + xHe + xNH3 )
    params.xd = xd
    """
    
    if use_vulcan == 0:
        params.atm_chemistry_arrays = {}
        for x in atm_chemistry:
            params.atm_chemistry_arrays['%s' % x] = [] 
            # initialized to replace the element corresponding to the current
            # integration step in the slope function
  
    if condensation == True:
        
        # Solve ODE               
        moist_w_cond = []                                             # Initialize the solution
        int_slope = integrator(slope,np.log(atm.ps),atm.temp[0],-1)  # Create the integrator instance. Negative increment to go from ps to ptop < ps
        int_slope.setParams(params)                                   # Update parameters used in the slope function dT/dlnP
        index = 0         # To count the number of iterations
        params.index2 = 0 # To count the number of iterations in slope()
        while int_slope.x > np.log(atm.ptop):                         # Stop at p = ptop

            if index > 1:
                T   = moist_w_cond[-1][1]
                lnP = moist_w_cond[-1][0]
            else:
                T   = atm.temp[0]
                lnP = np.log(atm.ps) 

            for molecule in atm_chemistry:  
                params.atm_chemistry_arrays[molecule].append(atm_chemistry[molecule])     
                # Check condensation for each molecule 
                if atm_chemistry[molecule] > 0.: # Tdew requires a non zero pressure
                    if T <= Tdew(molecule,params.atm_chemistry_arrays[molecule][-1]/numpy.exp(lnP)): 
                        params.atm_chemistry_arrays[molecule][-1] = esat(molecule,T)/numpy.exp(lnP) 

                    # print(molecule, T)
                    params.index2 += 1 

            moist_w_cond.append(int_slope.next())                     # Execute the Runge-Kutta integrator, fill array of tuples

            # print(moist_w_cond)
            
            #print(int_slope.x)
            index += 1

            # print(index, params.index2)
                        
        moist_w_cond = numpy.array(moist_w_cond)          # Convert to numpy array. lnP accessed through moist_w_cond[:,0]
                                                          #                         T   accessed through moist_w_cond[:,1]                   
        print(index, params.index2)                                                        
        print(params.atm_chemistry_arrays['H2O'])
        print(len(moist_w_cond[:,0]), len(params.atm_chemistry_arrays['H2O']))
        print(sum(1 for i in params.atm_chemistry_arrays['H2O'] if i < 0.999))  # Number of condensed levels
        print(sum(1 for i in params.atm_chemistry_arrays['H2O'] if i == 0.999)) # Number of non-condensed levels
        
    # Plot results
    
    # pressure array from the integral
    p_plot = numpy.exp(moist_w_cond[:,0]) # atm.p*1e-5 # bar
    
    Partial = False
    
    if Partial == False: # Condensation curves for a one-species atmosphere       
        TdewH2O = [ Tdew( 'H2O', pressure ) for pressure in p_plot ]
        TdewCO2 = [ Tdew( 'CO2', pressure ) for pressure in p_plot ]
        TdewCH4 = [ Tdew( 'CH4', pressure ) for pressure in p_plot ]
        TdewCO  = [ Tdew( 'CO',  pressure ) for pressure in p_plot ]
        TdewN2  = [ Tdew( 'N2',  pressure ) for pressure in p_plot ]
        TdewO2  = [ Tdew( 'O2',  pressure ) for pressure in p_plot ]
        TdewH2  = [ Tdew( 'H2',  pressure ) for pressure in p_plot ]
        TdewHe  = [ Tdew( 'He',  pressure ) for pressure in p_plot ]
        TdewNH3 = [ Tdew( 'NH3', pressure ) for pressure in p_plot ]
    
    else:                # Condensation curves for a mixed atmosphere (0 if species not present)
        if max(numpy.array(params.atm_chemistry_arrays['H2O'])) != 0.:
            TdewH2O = [ Tdew( 'H2O', pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['H2O'])/p_plot ]
        else:
            TdewH2O = numpy.zeros(len(p_plot))
        if max(numpy.array(params.atm_chemistry_arrays['CO2'])) != 0.:
            TdewCO2 = [ Tdew( 'CO2', pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['CO2'])/p_plot ]
        else:
            TdewCO2 = numpy.zeros(len(p_plot))
        if max(numpy.array(params.atm_chemistry_arrays['CH4'])) != 0.:
            TdewCH4 = [ Tdew( 'CH4', pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['CH4'])/p_plot ]
        else:
            TdewCH4 = numpy.zeros(len(p_plot))
        if max(numpy.array(params.atm_chemistry_arrays['CO'])) != 0.:
            TdewCO  = [ Tdew( 'CO',  pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['CO'])/p_plot ]
        else:
            TdewCO = numpy.zeros(len(p_plot))
        if max(numpy.array(params.atm_chemistry_arrays['N2'])) != 0.:
            TdewN2  = [ Tdew( 'N2',  pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['N2'])/p_plot ]
        else:
            TdewN2 = numpy.zeros(len(p_plot))
        if max(numpy.array(params.atm_chemistry_arrays['O2'])) != 0.:
            TdewO2  = [ Tdew( 'O2',  pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['O2'])/p_plot ]
        else:
            TdewO2 = numpy.zeros(len(p_plot))
        if  max(numpy.array(params.atm_chemistry_arrays['H2'])) != 0.:
            TdewH2  = [ Tdew( 'H2',  pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['H2'])/p_plot ]
        else:
            TdewH2 = numpy.zeros(len(p_plot))
        if max(numpy.array(params.atm_chemistry_arrays['He'])) != 0.:
            TdewHe  = [ Tdew( 'He',  pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['He'])/p_plot ]
        else:
            TdewHe = numpy.zeros(atm.nlev)
        if max(numpy.array(params.atm_chemistry_arrays['NH3'])) != 0.:
            TdewNH3 = [ Tdew( 'NH3', pressure ) for pressure in numpy.array(params.atm_chemistry_arrays['NH3'])/p_plot ]
        else:
            TdewNH3 = numpy.zeros(len(p_plot))
        
    sns.set_style("ticks")
    sns.despine()

    ls_adiabat = 2
    ls_ind = 0.8
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

    # plt.figure(1)
    ax1.semilogy(atm.ts*(p_plot/atm.ps)**(2./7.),p_plot,'r',lw=ls_adiabat,label=r'No condensation')
    
    if condensation == True:
        #print(moist_w_cond)
    
        ax1.semilogy(moist_w_cond[:,1],p_plot,'b',lw=ls_adiabat,label=r'With condensation')
        
    ax1.semilogy(TdewH2O,p_plot,label=r'H$_2$O', lw=ls_ind, ls="--")
    ax1.semilogy(TdewCO2,p_plot,label=r'CO$_2$', lw=ls_ind, ls="--")
    ax1.semilogy(TdewH2,p_plot, label=r'H$_2$',  lw=ls_ind, ls="--")
    ax1.semilogy(TdewCH4,p_plot,label=r'CH$_4$', lw=ls_ind, ls="--")
    ax1.semilogy(TdewCO,p_plot, label=r'CO',     lw=ls_ind, ls="--")
    ax1.semilogy(TdewN2,p_plot, label=r'N$_2$',  lw=ls_ind, ls="--")
    ax1.semilogy(TdewO2,p_plot, label=r'O$_2$',  lw=ls_ind, ls="--")
    # ax1.semilogy(TdewHe,p_plot, label=r'He', lw=ls_ind, ls="--")
    # ax1.semilogy(TdewNH3,p_plot,label=r'NH$_3$', lw=ls_ind, ls="--")
    
    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ (K)')
    ax1.set_ylabel(r'Total pressure $P$ (bar)')
    ax1.set_title('Individual moist adiabats')
    ax1.legend(ncol=1)
    #ax1.set_xlim([0,np.max(atm.temp)])
    ax1.set_ylim([np.max(p_plot),np.min(p_plot)])
    # plt.show()
    # plt.savefig('general_adiabat_TP.pdf', bbox_inches = 'tight')

    # plt.figure(2)
    
    
    xplot_H2O = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['H2O']),p_plot,label=r'H$_2$O')
    xplot_CO2 = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['CO2']),p_plot,label=r'CO$_2$')
    xplot_H2  = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['H2']),p_plot,label=r'H$_2$')
    xplot_CH4 = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['CH4']),p_plot,label=r'CH$_4$')
    xplot_CO  = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['CO']),p_plot,label=r'CO')
    xplot_N2  = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['N2']),p_plot,label=r'N$_2$')
    xplot_O2  = ax2.semilogy(numpy.array(params.atm_chemistry_arrays['O2']),p_plot,label=r'O$_2$')
    
    # Constant abundances w/o condensation as comparison
    
    if condensation == True:
    
        if use_vulcan == 0:
            const_abundances_H2O = np.ones(len(p_plot))*atm_chemistry["H2O"]
            const_abundances_CO2 = np.ones(len(p_plot))*atm_chemistry["CO2"]
            const_abundances_H2  = np.ones(len(p_plot))*atm_chemistry["H2"]
            const_abundances_CH4 = np.ones(len(p_plot))*atm_chemistry["CH4"]
            const_abundances_CO  = np.ones(len(p_plot))*atm_chemistry["CO"]
            const_abundances_N2  = np.ones(len(p_plot))*atm_chemistry["N2"]
            const_abundances_O2  = np.ones(len(p_plot))*atm_chemistry["O2"]
            const_abundances_He  = np.ones(len(p_plot))*atm_chemistry["He"]
            const_abundances_NH3 = np.ones(len(p_plot))*atm_chemistry["NH3"]
        ax2.semilogy(const_abundances_H2O, p_plot, ls="--", color=xplot_H2O[0].get_color() )
        ax2.semilogy(const_abundances_CO2, p_plot, ls="--", color=xplot_CO2[0].get_color() )
        ax2.semilogy(const_abundances_H2,  p_plot, ls="--", color=xplot_H2[0].get_color() )
        ax2.semilogy(const_abundances_CH4, p_plot, ls="--", color=xplot_CH4[0].get_color() )
        ax2.semilogy(const_abundances_CO,  p_plot, ls="--", color=xplot_CO[0].get_color() )
        ax2.semilogy(const_abundances_N2,  p_plot, ls="--", color=xplot_N2[0].get_color() )
        ax2.semilogy(const_abundances_O2,  p_plot, ls="--", color=xplot_O2[0].get_color() )

    # ax2.semilogy(xHe_array,p_plot,label=r'He')
    # ax2.semilogy(xNH3_array,p_plot,label=r'NH$_3$')
    ax2.invert_yaxis()
    ax2.set_xlabel(r'Mixing ratio $X_{\mathrm{vol}}/X_{\mathrm{tot}}$')
    ax2.set_ylabel(r'Total pressure $P$ (bar)')
    ax2.set_title('Relative abundances with condensation')
    ax2.legend(ncol=1)
    ax2.set_xlim(left=0, right=1.0)
    ax2.set_ylim([np.max(p_plot),np.min(p_plot)])
    # plt.show()
    plt.savefig('./output/general_adiabat.pdf', bbox_inches = 'tight')
    #plt.close(fig)

    return moist_w_cond   

# Define init parameters if called standalone
atm_chemistry  = { 
                "H2O" : 0.999, 
                "NH3" : 0.0,
                "CO2" : 0.0, 
                "CH4" : 0.0,
                "O2"  : 0.0,
                "CO"  : 0.0,
                "N2"  : 0.0, 
                "H2"  : 0.0,                        
                "He"  : 0.0  
                }
atm            = atmos()
atm.ts         = 300          # K
atm.ps         = 1e+6          # Pa
atm.ptop       = atm.ps*1e-5   # Pa
set_pressure_array(atm)
atm.temp        = atm.ts*(atm.p/atm.p[0])**atm.Rcp
use_vulcan     = 0
condensation   = True
moist_w_cond = solve_general_adiabat(atm, atm_chemistry, use_vulcan, condensation)

