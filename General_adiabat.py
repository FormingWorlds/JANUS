# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 13:17:05 2019

@author: Ryan
"""

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

def esat(switch,T): # Saturation vapor pressure for gases considered, in Pa
    #assuming ideal gas law and constant latent heat
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
    
#Dew point temperature
def Tdew(switch,p):
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

# Molar latent heat (J mol-1) for gases considered in general_moist_adiabat()
def L_heat(switch,T): 
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
    
def cpv(switch): # Molar heat capacities for gases considered, in J.K-1.mol-1
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

def general_moist_adiabat(lnP,T,xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3):

    # Kludge for H2 esat>7200 for T>13.95K
    if T < 13.95: 
        T = 13.95
    
    # Molar heat capacity of the dry species
    cpd = cpv('N2') #1000.*phys.air.MolecularWeight*1.e-3 

    # Avoid devision by zero if atmosphere consists of condensibles only
    xd = np.max([xd,1e-3])
                         
    xv_cpv = xH2O * cpv('H2O') + xCO2 * cpv('CO2') + \
                                 xCH4 * cpv('CH4') + \
                                 xCO  * cpv('CO')  + \
                                 xN2  * cpv('N2')  + \
                                 xO2  * cpv('O2')  + \
                                 xH2  * cpv('H2')  + \
                                 xHe  * cpv('He')  + \
                                 xNH3 * cpv('NH3')
  
    first_term = (xH2O/xd) * ( L_heat('H2O',T) / (R_universal*T) )**2 + \
                 (xCO2/xd) * ( L_heat('CO2',T) / (R_universal*T) )**2 + \
                 (xCH4/xd) * ( L_heat('CH4',T) / (R_universal*T) )**2 + \
                 (xCO/xd)  * ( L_heat('CO',T)  / (R_universal*T) )**2 + \
                 (xN2/xd)  * ( L_heat('N2',T)  / (R_universal*T) )**2 + \
                 (xO2/xd)  * ( L_heat('O2',T)  / (R_universal*T) )**2 + \
                 (xH2/xd)  * ( L_heat('H2',T)  / (R_universal*T) )**2 + \
                 (xHe/xd)  * ( L_heat('He',T)  / (R_universal*T) )**2 + \
                 (xNH3/xd) * ( L_heat('NH3',T) / (R_universal*T) )**2
    
    quadr_term=(
                 (xH2O/xd) * L_heat('H2O',T) / (R_universal*T) + \
                 (xCO2/xd) * L_heat('CO2',T) / (R_universal*T) + \
                 (xCH4/xd) * L_heat('CH4',T) / (R_universal*T) + \
                 (xCO/xd)  * L_heat('CO',T)  / (R_universal*T) + \
                 (xN2/xd)  * L_heat('N2',T)  / (R_universal*T) + \
                 (xO2/xd)  * L_heat('O2',T)  / (R_universal*T) + \
                 (xH2/xd)  * L_heat('H2',T)  / (R_universal*T) + \
                 (xHe/xd)  * L_heat('He',T)  / (R_universal*T) + \
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
              (xCO/xd)  * L_heat('CO',T)  / (R_universal*T) + \
              (xN2/xd)  * L_heat('N2',T)  / (R_universal*T) + \
              (xO2/xd)  * L_heat('O2',T)  / (R_universal*T) + \
              (xH2/xd)  * L_heat('H2',T)  / (R_universal*T) + \
              (xHe/xd)  * L_heat('He',T)  / (R_universal*T) + \
              (xNH3/xd) * L_heat('NH3',T) / (R_universal*T)
        
    # Factor T to get dT/dlnP
    num    = (1.+ num_sum)*T 
    denom  = (1./R_universal)*(xd*cpd+xv_cpv)/(xd+sum_abundances) + \
             (first_term+quadr_term)/(1.+sum_ratio_abundances)
    
    dTdlnP = num/denom    

    return dTdlnP    

# the new solve_ivp fct doesn't have arg() yet like odeint, but it will come in version 1.4.

def set_pressure_array2(atm):

    rat             = (atm.ptop/atm.ps)**(1./atm.nlev)
    logLevels       = [atm.ps*rat**i for i in range(atm.nlev+1)]
    logLevels.reverse()
    levels          = [atm.ptop + i*(atm.ps-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl          = np.array(logLevels)
    atm.p           = (atm.pl[1:] + atm.pl[:-1]) / 2

    return atm

def set_pressure_array(atm):

    atm.p    = np.exp(np.linspace(math.log(atm.ps),math.log(atm.ptop),atm.nlev))

    return atm

def solve_general_adiabat(atm, atm_chemistry, use_vulcan):

    # Set up pressure array (a global)
    atm  = set_pressure_array(atm)

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

    if use_vulcan == 0:
        xH2O_array = np.ones(len(atm.p))*xH2O
        xCO2_array = np.ones(len(atm.p))*xCO2
        xCH4_array = np.ones(len(atm.p))*xCH4
        xCO_array  = np.ones(len(atm.p))*xCO
        xN2_array  = np.ones(len(atm.p))*xN2
        xO2_array  = np.ones(len(atm.p))*xO2
        xH2_array  = np.ones(len(atm.p))*xH2
        xHe_array  = np.ones(len(atm.p))*xHe
        xNH3_array = np.ones(len(atm.p))*xNH3
        xd_array   = np.ones(len(atm.p))*xd

    # Solve ODE with dry, constant abundances
    atm.temp = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True) 

    # Check condensation for each species, adjust xv, solve ODE with new abundances
    for i in range(len(atm.p)):

        # If water is present
        if xH2O_array[i] != 0.:

            # If condensation occurs
            if atm.temp[i] <= Tdew('H2O',atm.p[i]): 

                # At surface 
                if i == 0: 

                    # Replace dry cst abundance with saturation vapor pressure at surface
                    atm.ps    = atm.ps - xH2O*atm.ps + esat('H2O',atm.ts)
                    
                    # Update the pressure array with the new ps. Keep ptop cst?
                    atm       = set_pressure_array(atm)

                # Update the new abundances
                xH2O          = esat('H2O',atm.temp[i])/atm.p[i] 
                xd            = 1. - (xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xH2O_array[i] = xH2O      
                xd_array[i]   = xd

                # Recalculate the moist adiabat with new abundances
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)

        if xNH3_array[i] != 0.:
            if atm.temp[i] < Tdew('NH3',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xNH3*atm.ps + esat('NH3',atm.ts)
                    atm       = set_pressure_array(atm)
                xNH3          = esat('NH3',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)       
                xNH3_array[i] = xNH3
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
           
        if xCO2_array[i] != 0.:
            if atm.temp[i] < Tdew('CO2',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xCO2*atm.ps + esat('CO2',atm.ts)
                    atm       = set_pressure_array(atm)
                xCO2          = esat('CO2',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xCO2_array[i] = xCO2
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)

        if xCH4_array[i] != 0.:
            if atm.temp[i] < Tdew('CH4',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xCH4*atm.ps + esat('CH4',atm.ts)
                    atm       = set_pressure_array(atm)
                xCH4          = esat('CH4',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xCH4_array[i] = xCH4
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)

        if xO2_array[i] != 0.:
            if atm.temp[i] < Tdew('O2',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xO2*atm.ps + esat('O2',atm.ts)
                    atm       = set_pressure_array(atm)
                xO2           = esat('O2',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xO2_array[i]  = xO2
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
            
        if xCO_array[i] != 0.:
            if atm.temp[i] < Tdew('CO',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xCO*atm.ps + esat('CO',atm.ts)
                    atm       = set_pressure_array(atm)
                xCO           = esat('CO',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xCO_array[i]  = xCO
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)

        if xN2_array[i] != 0.:
            if atm.temp[i] < Tdew('N2',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xN2*atm.ps + esat('N2',atm.ts)
                    atm       = set_pressure_array(atm)
                xN2           = esat('N2',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xN2_array[i]  = xN2
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
            
                # Kludge for H2 esat > 7200 for T > 13.95K
                if atm.temp[i] < 13.95: 
                    atm.temp[i] = 13.95
             
        if xH2_array[i] != 0.:
            if atm.temp[i] < Tdew('H2',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xH2*atm.ps + esat('H2',atm.ts)
                    atm       = set_pressure_array(atm)
                xH2           = esat('H2',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
                xH2_array[i]  = xH2
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
            
        if xHe_array[i] != 0.:
            if atm.temp[i] < Tdew('He',atm.p[i]):
                if i == 0:
                    atm.ps    = atm.ps - xHe*atm.ps + esat('He',atm.ts)
                    atm       = set_pressure_array(atm)
                xHe           = esat('He',atm.temp[i])/atm.p[i]
                xd            = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)       
                xHe_array[i]  = xHe
                xd_array[i]   = xd
                atm.temp      = odeint(general_moist_adiabat,atm.ts,np.log(atm.p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
            
    # Plot results
            
    TdewH2O = [ Tdew( 'H2O', pressure ) for pressure in atm.p ]
    TdewCO2 = [ Tdew( 'CO2', pressure ) for pressure in atm.p ]
    TdewCH4 = [ Tdew( 'CH4', pressure ) for pressure in atm.p ]
    TdewCO  = [ Tdew( 'CO',  pressure ) for pressure in atm.p ]
    TdewN2  = [ Tdew( 'N2',  pressure ) for pressure in atm.p ]
    TdewO2  = [ Tdew( 'O2',  pressure ) for pressure in atm.p ]
    TdewH2  = [ Tdew( 'H2',  pressure ) for pressure in atm.p ]
    TdewHe  = [ Tdew( 'He',  pressure ) for pressure in atm.p ]
    TdewNH3 = [ Tdew( 'NH3', pressure ) for pressure in atm.p ]

    sns.set_style("ticks")
    sns.despine()

    p_plot = atm.p*1e-5 # bar

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

    # plt.figure(1)
    ax1.semilogy(atm.temp,p_plot,'r',linewidth=3,label=r'Combined')
    ax1.semilogy(TdewH2O,p_plot,label=r'H$_2$O',ls="--")
    ax1.semilogy(TdewCO2,p_plot,label=r'CO$_2$',ls="--")
    ax1.semilogy(TdewH2,p_plot,label=r'H$_2$',ls="--")
    ax1.semilogy(TdewCH4,p_plot,label=r'CH$_4$',ls="--")
    ax1.semilogy(TdewCO,p_plot,label=r'CO',ls="--")
    ax1.semilogy(TdewN2,p_plot,label=r'N$_2$',ls="--")
    ax1.semilogy(TdewO2,p_plot,label=r'O$_2$',ls="--")
    # ax1.semilogy(TdewHe,p_plot,label=r'He')
    # ax1.semilogy(TdewNH3,p_plot,label=r'NH$_3$')
    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ (K)')
    ax1.set_ylabel(r'Total pressure $P$ (bar)')
    ax1.set_title('Individual moist adiabats')
    ax1.legend(ncol=1)
    ax1.set_xlim([0,np.max(atm.temp)])
    ax1.set_ylim([np.max(p_plot),np.min(p_plot)])
    # plt.show()
    # plt.savefig('general_adiabat_TP.pdf', bbox_inches = 'tight')

    # plt.figure(2)
    ax2.semilogy(xH2O_array,p_plot,label=r'H$_2$O')
    ax2.semilogy(xCO2_array,p_plot,label=r'CO$_2$')
    ax2.semilogy(xH2_array,p_plot,label=r'H$_2$')
    ax2.semilogy(xCH4_array,p_plot,label=r'CH$_4$')
    ax2.semilogy(xCO_array,p_plot,label=r'CO')
    ax2.semilogy(xN2_array,p_plot,label=r'N$_2$')
    ax2.semilogy(xO2_array,p_plot,label=r'O$_2$') 
    # ax2.semilogy(xHe_array,p_plot,label=r'He')
    # ax2.semilogy(xNH3_array,p_plot,label=r'NH$_3$')
    ax2.invert_yaxis()
    ax2.set_xlabel(r'Mixing ratio $X_{\mathrm{vol}}/X_{\mathrm{tot}}$')
    ax2.set_ylabel(r'Total pressure $P$ (bar)')
    ax2.set_title('Relative abundances')
    ax2.legend(ncol=1)
    ax2.set_xlim(left=0)
    ax2.set_ylim([np.max(p_plot),np.min(p_plot)])
    # plt.show()
    plt.savefig('./output/general_adiabat.pdf', bbox_inches = 'tight')
    plt.close(fig)

    """
    plt.semilogy(Moist_adiabat,p)
    plt.semilogy(TdewH2O,p,label='H2O')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    """

    return atm.temp    

# Define init parameters if called standalone
atm_chemistry  = { 
                "H2O" : 0.5, 
                "CO2" : 0.3, 
                "H2"  : 0.0, 
                "CH4" : 0.0,
                "CO"  : 0.0,  
                "N2"  : 0.0, 
                "O2"  : 0.0,
                "He"  : 0.0,  
                "NH3" : 0.0
                }
atm            = atmos()
atm.ts         = 600          # K
atm.ps         = 1e+6          # Pa
atm.ptop       = atm.ps*1e-5   # Pa
use_vulcan     = 0
atm_moist_temp = solve_general_adiabat(atm, atm_chemistry, use_vulcan)

