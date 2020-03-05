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

def esat(switch,T): 
    """Saturation vapor pressure [Pa] given a temperature T [K]. Assuming the ideal gas law and a constant latent heat. Select the molecule of interest with the switch argument (a string)."""
    if switch == 'H2O':
        e=phys.satvps_function(phys.water)
    if switch == 'CH4':
        e=phys.satvps_function(phys.methane)
    if switch == 'CO2':
        e=phys.satvps_function(phys.co2)
    if switch == 'CO':
        e=phys.satvps_function(phys.co)
    if switch == 'N2':
        e=phys.satvps_function(phys.n2)
    if switch == 'O2':
        e=phys.satvps_function(phys.o2)
    if switch == 'H2':
        e=phys.satvps_function(phys.h2)
    if switch == 'He':
        e=phys.satvps_function(phys.he)
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

    print(switch, T, T_triple, T_crit, L_heat)

    return L_heat  
    
## Dew point temperature [K] array for given pressure array [Pa] and surface T [K]. 
## Select the molecule of interest with the switch argument (a string)
def Tdew(switch, prs, T_surf, L_use): 

    # Calculate dew-point for each pressure
    T_dew = T_surf / ( 1. - ( T_surf * R_universal / L_use ) * np.log( prs / np.max(prs) ) )

    # Re-calc with correct L_heat
    for idx, T in enumerate(T_dew):
        L_use = L_heat( switch, T )
        T_dew[idx] = T_surf / ( 1. - ( T_surf * R_universal / L_use ) * np.log( prs[idx] / np.max(prs) ) )

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


def set_pressure_array(atm):
    rat       = (atm.ptop/atm.ps)**(1./atm.nlev)
    logLevels = [atm.ps*rat**i for i in range(atm.nlev+1)]
    levels    = [atm.ptop + i*(atm.ps-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl    = np.array(logLevels)
    atm.p     = (atm.pl[1:] + atm.pl[:-1]) / 2
    return atm
                        
def solve_general_adiabat(atm, atm_chemistry, use_vulcan, condensation):
    """Builds the generalized moist adiabat from slope(lnP,T,params) and plots the adiabat along with the relative abundances."""   

    sns.set_style("ticks")
    sns.despine()

    ls_adiabat = 2
    ls_ind = 0.8

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))

    for mode in [ "original" ]: # "hybrid", , "ray1"

        params = Dummy() # initialize parameter object  
    
        if use_vulcan == 0:
            params.atm_chemistry_arrays = {}
            for x in atm_chemistry:
                params.atm_chemistry_arrays['%s' % x] = [] 
              
        # --------------------------------- Initialisations ---------------------------------          
        moist_w_cond = [] #[tuple([np.log(atm.ps), atm.ts])]     # Initialize the tuple solution
        pL = []                                                  # Initialize the final pressure array
        psatL = []
        Rcp = []                                                 # Initialize the exponent of the dry adiabat
        TL = []                                                  # Initialize the final temperature array
        step = -.1                                               # Negative increment to go from ps to ptop < ps
        int_slope = integrator(slope,np.log(atm.ps),np.log(atm.ts),step) # Create the integrator instance.
        int_slopeRay = integrator(slopeRay,np.log(atm.ps),np.log(atm.ts),step)
        int_slope.setParams(params)                              # Update parameters used in the slope function dT/dlnP
        index = 0    


        logT = int_slope.y
        logP = int_slope.x
        Temperature = np.exp(logT)
        Pressure = np.exp(logP)


                                                    # To count the number of iterations
        
        # ----------------------------------- Integration -----------------------------------
        """ Original attempt """
        if mode == "original":
            while Pressure > atm.ptop:                    # Start at ln(ps), stop at ln(ptop)
                
                for molecule in atm_chemistry:                       # Loop on the molecules 
                    params.atm_chemistry_arrays[molecule].append(atm_chemistry[molecule]) # Abundance at p[0]                 
                    if atm_chemistry[molecule] > 0.:                                      # Tdew requires a non zero pressure
                        p_molecule = params.atm_chemistry_arrays[molecule][-1]*Pressure # partial pressure
                        p_sat = esat(molecule,Temperature)
                        #print(molecule,p_molecule,p_sat,params.atm_chemistry_arrays[molecule][-1])
                        #print(params.atm_chemistry_arrays[molecule])
                        # if Temperature <= Tdew(molecule,Pressure,atm.ts): # Check condensation for each molecule
                        
                        # If int_slope.y < Tdew, how can numpy.exp(int_slope.x) < p_sat?
                        params.atm_chemistry_arrays[molecule][-1] = np.min([p_molecule,p_sat])/Pressure # min to avoid p_sat > p_molecule       
                        # params.atm_chemistry_arrays[molecule][-1] = p_sat/Pressure
                            
                moist_w_cond.append(int_slope.next()) # Execute the Runge-Kutta integrator, fill array of tuples
                pL.append(numpy.exp(int_slope.x))
                TL.append(numpy.exp(int_slope.y))

                logT = int_slope.y
                logP = int_slope.x
                Temperature = np.exp(logT)
                Pressure = np.exp(logP)
                Rcp.append(R_universal/params.general_cp)


                index += 1

            # print(atm.ps, atm.ts)
            # print(esat('H2O',atm.ts), np.max(params.atm_chemistry_arrays['H2O'])*np.max(pL))
            # print(np.max(pL), np.max(TL))

        
        
        # """ Basically the same thing, but p_total = p_air + p_sat 
        if mode == "hybrid":
            # Difference between Ray's slope and ours: abundances taken from last element of atm_chemistry_arrays, 
            # so they need at least one element before first integration step.
            for molecule in atm_chemistry:
                    params.atm_chemistry_arrays[molecule].append(atm_chemistry[molecule])
                    #print(params.atm_chemistry_arrays)
            while int_slope.x > np.log(atm.ptop):        
                index += 1
                ans = int_slope.next()
                moist_w_cond.append(ans)
                pa = math.exp(ans[0])
                T = math.exp(ans[1])
                p = pa+esat('H2O',T)
                pL.append(p)
                #print(index)
                #print(esat('H2O',T))
                for molecule in atm_chemistry:
                    #if atm_chemistry[molecule] > 0.:  
                    #if numpy.exp(int_slope.y) <= Tdew(molecule,numpy.exp(int_slope.x)):
                    if index == 1 :
                    #print(index)
                        if atm_chemistry[molecule] > 0.:
                            params.atm_chemistry_arrays[molecule][0] = esat(molecule,T)/p
                        else:
                            params.atm_chemistry_arrays[molecule][0] = atm_chemistry[molecule]
                        #print(params.atm_chemistry_arrays)
                    else:
                        if atm_chemistry[molecule] > 0.:                  
                            #print(index)
                            params.atm_chemistry_arrays[molecule].append(esat(molecule,T)/p)
                            #print(params.atm_chemistry_arrays)
                            #params.atm_chemistry_arrays[molecule][-1] = esat(molecule,numpy.exp(int_slope.y))/numpy.exp(int_slope.x)
                        else:
                            params.atm_chemistry_arrays[molecule].append(atm_chemistry[molecule])
                    
                TL.append(T) 
        
        #print(params.atm_chemistry_arrays['H2O'])
        #print(pL)
        # """


        # """ Exact same, but with Ray's slope
        if mode == "ray1":
            while int_slopeRay.x > np.log(atm.ptop):                 
                for molecule in atm_chemistry:
                    params.atm_chemistry_arrays[molecule].append(atm_chemistry[molecule]) 
                ans = int_slopeRay.next()
                moist_w_cond.append(ans)
                pa = math.exp(ans[0])
                T = math.exp(ans[1])
                p = pa+esat('H2O',T)
                pL.append(p)
                psatL.append(esat('H2O',T))
                for molecule in atm_chemistry:
                    if atm_chemistry[molecule] > 0.:    
                        #params.atm_chemistry_arrays[molecule].append(esat('H2O',T)/p)
                        params.atm_chemistry_arrays[molecule][-1] = esat(molecule,T)/p            
                TL.append(T) 
                index += 1

            # print(len(params.atm_chemistry_arrays['H2O']))
        # """
        
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
        # print(len(pL), len(TL), len(abundance_arrayH2O))
        # print(abundance_arrayH2O)
        
        # print(sum(1 for i in params.atm_chemistry_arrays['H2O'] if i < 0.999))  # Number of condensed levels
        # print(sum(1 for i in params.atm_chemistry_arrays['H2O'] if i == 0.999)) # Number of non-condensed levels
            
        # ------------------------------------ Plotting -------------------------------------

        # # pressure array from the integral
        # p_int = numpy.exp(moist_w_cond[:,0]) # atm.p*1e-5 # bar
        # # pressure array for plotting

        # Pressure array for interpolation scheme
        # p_plot = np.logspace(-5,5,100)
        p_plot = np.exp(np.linspace(np.log(atm.ptop),np.log(atm.ps),100))

        # Execute Ray's function
        if mode == "original":
            # Compute Ray's moist adiabat (valid for a single condensible gas)
            moist_ray  = phys.MoistAdiabat(phys.H2O,phys.N2)
            p_noncondensible = atm.ps*(1.-atm_chemistry["H2O"])
            # print(atm.ts, esat('H2O',atm.ts))
            # print("p_noncondensible: ", p_noncondensible)
            p_ray,T_ray,molarCon,massCon = moist_ray(p_noncondensible,atm.ts,np.min(p_plot))
            p_ray_interp,T_ray_interp,molarCon_interp,massCon_interp = moist_ray(p_noncondensible,atm.ts,np.min(p_plot),p_plot)
            # print(p_ray)
            # print(p_plot)
            #print(molarCon_interp)
            # p_ray_interp is a copy of p_plot
            # print(p_ray_interp)
            # print(T_ray_interp)
            # print(molarCon_interp)
            # print(massCon_interp)
            # print(p_noncondensible)
            #print([esat('H2O',T_ray_interp[i])/p_ray_interp[i] for i in range(len(p_ray_interp))])

            # Plot Ray's H2O moist adiabat function
            ax1.semilogy(T_ray_interp,p_plot,lw=ls_adiabat, color="gray", ls="--",label=r'p$_{non-cond.}$ = '+"{:.2f}".format(p_noncondensible)+' Pa')
            ax2.semilogy(molarCon_interp,p_plot,lw=ls_adiabat, color="gray", ls="--",label=r"Ray's function")
        
        
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
            
            T_interp                  = numpy.array([T1(pp) for pp in pL])
            abundance_arrayH2O_interp = numpy.array([abundance_arrayH2O1(pp) for pp in pL])
            abundance_arrayNH3_interp = numpy.array([abundance_arrayNH31(pp) for pp in pL])
            abundance_arrayCO2_interp = numpy.array([abundance_arrayCO21(pp) for pp in pL])
            abundance_arrayCH4_interp = numpy.array([abundance_arrayCH41(pp) for pp in pL])
            abundance_arrayCO_interp  = numpy.array([abundance_arrayCO1(pp) for pp in pL])
            abundance_arrayN2_interp  = numpy.array([abundance_arrayN21(pp) for pp in pL])
            abundance_arrayO2_interp  = numpy.array([abundance_arrayO21(pp) for pp in pL])
            abundance_arrayH2_interp  = numpy.array([abundance_arrayH21(pp) for pp in pL])
            abundance_arrayHe_interp  = numpy.array([abundance_arrayHe1(pp) for pp in pL])

            pL = p_plot


        # Plot single-species dew points
        for vol in [ 'H2O' ]:
            L_use = L_heat( vol, atm.ts )
            Tdew_array = Tdew( vol, pL, atm.ts, L_use )
            ax1.semilogy( Tdew_array, pL, label=vol_latex[vol]+' dew-point', lw=ls_ind, ls=":", color=vol_colors[vol+"_1"])

        #ax1.semilogy(moist_w_cond[:,1],p_plot,color="red",lw=ls_adiabat,label=r'Moist adiabat')

        if mode ==  "original":
            color = "red"
        if mode ==  "hybrid":
            color = "blue"
        if mode ==  "ray1":
            color = "green"
        label = mode
            
            
            #ax1.semilogy(T_ray,p_ray,lw=ls_adiabat, color="blue", ls="-")

        if mode == "original":

            ax1.semilogy(atm.ts*(pL/atm.ps)**Rcp,pL,color='black', ls="--",lw=ls_adiabat,label=r'Dry adiabat')
            # print(Rcp)
            
        alpha=0.0
        
        if Interpolate:
            xplot_H2O = ax2.semilogy(abundance_arrayH2O_interp,pL,label=mode+r': H$_2$O', color=color)
            xplot_CO2 = ax2.semilogy(abundance_arrayCO2_interp,pL, alpha=alpha)#,label=r'CO$_2$', alpha=alpha)
            xplot_H2  = ax2.semilogy(abundance_arrayH2_interp,pL, alpha=alpha)#,label=r'H$_2$', alpha=alpha)
            xplot_CH4 = ax2.semilogy(abundance_arrayCH4_interp,pL, alpha=alpha)#,label=r'CH$_4$', alpha=alpha)
            xplot_CO  = ax2.semilogy(abundance_arrayCO_interp,pL, alpha=alpha)#,label=r'CO', alpha=alpha)
            xplot_N2  = ax2.semilogy(abundance_arrayN2_interp,pL, alpha=alpha)#,label=r'N$_2$', alpha=alpha)
            xplot_O2  = ax2.semilogy(abundance_arrayO2_interp,pL, alpha=alpha)#,label=r'O$_2$', alpha=alpha)
        else:
            xplot_H2O = ax2.semilogy(abundance_arrayH2O,pL,label=mode+r': H$_2$O', color=color)
            xplot_CO2 = ax2.semilogy(abundance_arrayCO2,pL, alpha=alpha)#,label=r'CO$_2$', alpha=alpha)#
            xplot_H2  = ax2.semilogy(abundance_arrayH2,pL, alpha=alpha)#,label=r'H$_2$', alpha=alpha)
            xplot_CH4 = ax2.semilogy(abundance_arrayCH4,pL, alpha=alpha)#,label=r'CH$_4$', alpha=alpha)
            xplot_CO  = ax2.semilogy(abundance_arrayCO,pL, alpha=alpha)#,label=r'CO', alpha=alpha)
            xplot_N2  = ax2.semilogy(abundance_arrayN2,pL, alpha=alpha)#,label=r'N$_2$', alpha=alpha)
            xplot_O2  = ax2.semilogy(abundance_arrayO2,pL, alpha=alpha)#,label=r'O$_2$', alpha=alpha)

        # Calculated moist adiabat
        if Interpolate:
            ax1.semilogy(T_interp,pL,color=color,lw=ls_adiabat,label="Moist adiabat",alpha=0.) # interpolated
        else:
            ax1.semilogy(TL,pL,color=color,lw=ls_adiabat,label="Moist adiabat",alpha=0.) # non-interpolated

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

    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ (K)')
    ax1.set_ylabel(r'Total pressure $P$ (Pa)')
    ax1.set_title('Adiabats + individual Clausius-Clapeyron slopes')
    ax1.legend(ncol=1)
    # ax1.set_xlim([0,np.max(atm.ts)])
    # ax1.set_ylim([1e5,1e-5])
    # ax1.set_ylim([np.max(p_plot),np.min(p_plot)])

    # ax2.semilogy(xHe_array,p_plot,label=r'He')
    # ax2.semilogy(xNH3_array,p_plot,label=r'NH$_3$')
    ax2.invert_yaxis()
    ax2.set_title('Mixing ratios')
    ax2.set_xlabel(r'Molar concentration $X_{\mathrm{cond.}}/X_{total}$')
    ax2.set_ylabel(r'Total pressure $P$ (Pa)')
    # ax2.set_title('Relative abundances with condensation')
    ax2.legend(ncol=1)
    ax2.set_xlim(left=-0.05, right=1.05)
    ax2.set_ylim([1e5,1e-5])
    
    # plt.show()
    plt.savefig('./output/general_adiabat.pdf', bbox_inches = 'tight')
    #plt.close(fig)

    return moist_w_cond   

# Define init parameters if called standalone
atm_chemistry  = { 
                "H2O" : 0.5, 
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
atm.ts         = 500.           # K
atm.ps         = 1e+5           # Pa
atm.ptop       = atm.ps*1e-10   # Pa
set_pressure_array(atm)
atm.temp       = atm.ts*(atm.p/atm.p[0])**atm.Rcp
use_vulcan     = 0
condensation   = True
moist_w_cond = solve_general_adiabat(atm, atm_chemistry, use_vulcan, condensation)

