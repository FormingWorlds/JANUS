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

#--------- Importing thermodynamical properties of gases -----------

R_universal=8.31446261815324 # Universal gas constant, J.K-1.mol-1, should probably go elsewhere since it's a constant    

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
    
def L_heat(switch,T): # Molar latent heat for gases considered, in J.mol-1. Used in General_moist_adiabat()
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

def General_moist_adiabat(lnP,T,xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3):

    if T<13.95: # Kludge, for H2 esat>7200 for T>13.95K which seems to be a hard limit
        T=13.95
    
    cpd=cpv('N2') #1000.*phys.air.MolecularWeight*1.e-3 # molar heat capacity of the dry species
                         
    xv_cpv=xH2O*cpv('H2O') + xCO2*cpv('CO2') +\
                             xCH4*cpv('CH4') +\
                             xCO*cpv('CO')   +\
                             xN2*cpv('N2')   +\
                             xO2*cpv('O2')   +\
                             xH2*cpv('H2')   +\
                             xHe*cpv('He')   +\
                             xNH3*cpv('NH3')
  
    first_term=(xH2O/xd)*(L_heat('H2O',T)/(R_universal*T))**2 +\
               (xCO2/xd)*(L_heat('CO2',T)/(R_universal*T))**2 +\
               (xCH4/xd)*(L_heat('CH4',T)/(R_universal*T))**2 +\
               (xCO/xd)*(L_heat('CO',T)/(R_universal*T))**2   +\
               (xN2/xd)*(L_heat('N2',T)/(R_universal*T))**2   +\
               (xO2/xd)*(L_heat('O2',T)/(R_universal*T))**2   +\
               (xH2/xd)*(L_heat('H2',T)/(R_universal*T))**2   +\
               (xHe/xd)*(L_heat('He',T)/(R_universal*T))**2   +\
               (xNH3/xd)*(L_heat('NH3',T)/(R_universal*T))**2
    
    quadr_term=(
               (xH2O/xd)*L_heat('H2O',T)/(R_universal*T) +\
               (xCO2/xd)*L_heat('CO2',T)/(R_universal*T) +\
               (xCH4/xd)*L_heat('CH4',T)/(R_universal*T) +\
               (xCO/xd)*L_heat('CO',T)/(R_universal*T)   +\
               (xN2/xd)*L_heat('N2',T)/(R_universal*T)   +\
               (xO2/xd)*L_heat('O2',T)/(R_universal*T)   +\
               (xH2/xd)*L_heat('H2',T)/(R_universal*T)   +\
               (xHe/xd)*L_heat('He',T)/(R_universal*T)   +\
               (xNH3/xd)*L_heat('NH3',T)/(R_universal*T)             
                                                       )**2
        
    sum_abundances=xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3
                                             
    sum_ratio_abundances=xH2O/xd + xCO2/xd +\
                                   xCH4/xd +\
                                   xCO/xd  +\
                                   xN2/xd  +\
                                   xO2/xd  +\
                                   xH2/xd  +\
                                   xHe/xd  +\
                                   xNH3/xd
                                                      
    num_sum=(xH2O/xd)*L_heat('H2O',T)/(R_universal*T) +\
            (xCO2/xd)*L_heat('CO2',T)/(R_universal*T) +\
            (xCH4/xd)*L_heat('CH4',T)/(R_universal*T) +\
            (xCO/xd)*L_heat('CO',T)/(R_universal*T)   +\
            (xN2/xd)*L_heat('N2',T)/(R_universal*T)   +\
            (xO2/xd)*L_heat('O2',T)/(R_universal*T)   +\
            (xH2/xd)*L_heat('H2',T)/(R_universal*T)   +\
            (xHe/xd)*L_heat('He',T)/(R_universal*T)   +\
            (xNH3/xd)*L_heat('NH3',T)/(R_universal*T)
        
    num=(1.+ num_sum)*T # Factor T to get dT/dlnP
    denom=(1./R_universal)*(xd*cpd+xv_cpv)/(xd+sum_abundances) + (first_term+quadr_term)/(1.+sum_ratio_abundances)
    dTdlnP=num/denom       
    return dTdlnP    

# the new solve_ivp fct doesn't have arg() yet like odeint, but it will come in version 1.4.
                   
# Ground temperature
T0   = 373.15

# Initialising the relative abundances # Careful about condensing things that are not there...
xH2O = 0.2
xCO2 = 0.2
xCH4 = 0.
xCO  = 0.
xN2  = 0.
xO2  = 0.
xH2  = 0.
xHe  = 0.
xNH3 = 0.
xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)

#lnP points
nlev = 100         # number of pressure levels
ps   = 1e7         # Surface pressure in Pa
ptop = 1.e-5*ps    # Pressure at the TOA in Pa
p    = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev)) # Logarithmic pressure grid scale

# Initialising the arrays for plotting
xH2O_array = np.ones(len(p))*xH2O
xCO2_array = np.ones(len(p))*xCO2
xCH4_array = np.ones(len(p))*xCH4
xCO_array  = np.ones(len(p))*xCO
xN2_array  = np.ones(len(p))*xN2
xO2_array  = np.ones(len(p))*xO2
xH2_array  = np.ones(len(p))*xH2
xHe_array  = np.ones(len(p))*xHe
xNH3_array = np.ones(len(p))*xNH3
xd_array   = np.ones(len(p))*xd
 
#solve ODE with dry, constant abundances
Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True) 
#Moist_adiabat=solve_ivp(lambda lnP,T: General_moist_adiabat(lnP,T,xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),(np.log(ps),np.log(ptop)),T0*np.ones(len(p)))
# This works only if T is an array in General_moist_adiabat().

#Check condensation for each species, adjust xv if needed and solve ODE with the new abundances
for i in range(len(p)):
    if xH2O_array[i] != 0.: # If water is present
        if Moist_adiabat[i]<=Tdew('H2O',p[i]): # If condensation occurs
            #print("condensation")
            if i==0: # 0 is surface level
                ps = ps - xH2O*ps + esat('H2O',T0) # We replace the dry cst abundance with the saturation vapor pressure at the surface
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev)) # Update the pressure array with the new ps. keeping ptop cst?
            xH2O=esat('H2O',Moist_adiabat[i])/p[i] # Update the new abundance where needed
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3) # Update xd
            xH2O_array[i]=xH2O      
            xd_array[i]=xd
            #xH2O_array = np.where(xH2O_array>1.,1.,xH2O_array) # Kludge 
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array)) 
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True) 

    if xNH3_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('NH3',p[i]):
            if i==0:
                ps = ps - xNH3*ps + esat('NH3',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xNH3=esat('NH3',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)       
            xNH3_array[i]=xNH3
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
       
    if xCO2_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('CO2',p[i]):
            if i==0:
                ps = ps - xCO2*ps + esat('CO2',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xCO2=esat('CO2',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
            xCO2_array[i]=xCO2
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True) 

    if xCH4_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('CH4',p[i]):
            if i==0:
                ps = ps - xCH4*ps + esat('CH4',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xCH4=esat('CH4',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
            xCH4_array[i]=xCH4
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True) 

    if xO2_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('O2',p[i]):
            if i==0:
                ps = ps - xO2*ps + esat('O2',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xO2=esat('O2',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
            xO2_array[i]=xO2
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
        
    if xCO_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('CO',p[i]):
            if i==0:
                ps = ps - xCO*ps + esat('CO',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xCO=esat('CO',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
            xCO_array[i]=xCO
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True) 

    if xN2_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('N2',p[i]):
            if i==0:
                ps = ps - xN2*ps + esat('N2',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xN2=esat('N2',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
            xN2_array[i]=xN2
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
        
            if Moist_adiabat[i]<13.95: # Kludge, for H2 esat>7200 for T>13.95K which seems to be a hard limit
                Moist_adiabat[i]=13.95
         
    if xH2_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('H2',p[i]):
            if i==0:
                ps = ps - xH2*ps + esat('H2',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xH2=esat('H2',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)
            xH2_array[i]=xH2
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
        
    if xHe_array[i] != 0.:
        if Moist_adiabat[i]<Tdew('He',p[i]):
            if i==0:
                ps = ps - xHe*ps + esat('He',T0)
                p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev))
            xHe=esat('He',Moist_adiabat[i])/p[i]
            xd   = 1.-(xH2O+xCO2+xCH4+xCO+xN2+xO2+xH2+xHe+xNH3)       
            xHe_array[i]=xHe
            xd_array[i]=xd
            #Moist_adiabat[i]=odeint(General_moist_adiabat,T0,np.log(p[i]),args=(xd_array[i],xH2O_array[i],xCO2_array[i],xCH4_array[i],xCO_array[i],xN2_array[i],xO2_array[i],xH2_array[i],xHe_array[i],xNH3_array[i]))
            #Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd_array,xH2O_array,xCO2_array,xCH4_array,xCO_array,xN2_array,xO2_array,xH2_array,xHe_array,xNH3_array))
            Moist_adiabat=odeint(General_moist_adiabat,T0,np.log(p),args=(xd,xH2O,xCO2,xCH4,xCO,xN2,xO2,xH2,xHe,xNH3),tfirst=True)
        
        
#plot results
        
TdewH2O = [ Tdew('H2O',pp) for pp in p ]
TdewCO2 = [ Tdew('CO2',pp) for pp in p ]
TdewCH4 = [ Tdew('CH4',pp) for pp in p ]
TdewCO  = [ Tdew('CO',pp) for pp in p ]
TdewN2  = [ Tdew('N2',pp) for pp in p ]
TdewO2  = [ Tdew('O2',pp) for pp in p ]
TdewH2  = [ Tdew('H2',pp) for pp in p ]
TdewHe  = [ Tdew('He',pp) for pp in p ]
TdewNH3 = [ Tdew('NH3',pp) for pp in p ]

plt.figure(1)
plt.semilogy(Moist_adiabat,p,'r',linewidth=4)
plt.semilogy(TdewH2O,p,label='H2O')
plt.semilogy(TdewCO2,p,label='CO2')
plt.semilogy(TdewCH4,p,label='CH4')
plt.semilogy(TdewCO,p,label='CO')
plt.semilogy(TdewN2,p,label='N2')
plt.semilogy(TdewO2,p,label='O2')
plt.semilogy(TdewH2,p,label='H2')
plt.semilogy(TdewHe,p,label='He')
plt.semilogy(TdewNH3,p,label='NH3')
plt.gca().invert_yaxis()
plt.xlabel('T (K)')
plt.ylabel('P (Pa)')
plt.title('Individual moist adiabats')
plt.legend()
plt.show()
plt.savefig('general_adiabat_TP.pdf', bbox_inches = 'tight')

plt.figure(2)
plt.semilogy(xH2O_array,p,label='H2O relative abundance')
plt.semilogy(xCO2_array,p,label='CO2 relative abundance')
plt.semilogy(xCH4_array,p,label='CH4 relative abundance')
plt.semilogy(xCO_array,p,label='CO relative abundance')
plt.semilogy(xN2_array,p,label='N2 relative abundance')
plt.semilogy(xO2_array,p,label='O2 relative abundance')
plt.semilogy(xH2_array,p,label='H2 relative abundance')
plt.semilogy(xHe_array,p,label='He relative abundance')
plt.semilogy(xNH3_array,p,label='NH3 relative abundance')
plt.gca().invert_yaxis()
plt.xlabel('xv')
plt.ylabel('P (Pa)')
plt.title('Relative abundances')
plt.legend()
plt.show()
plt.savefig('general_adiabat_XP.pdf', bbox_inches = 'tight')

"""
plt.semilogy(Moist_adiabat,p)
plt.semilogy(TdewH2O,p,label='H2O')
plt.gca().invert_yaxis()
plt.legend()
plt.show()
"""