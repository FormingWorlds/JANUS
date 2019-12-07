# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 13:17:05 2019

@author: Ryan
"""

import numpy as np
import math,phys
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#--------- Importing thermodynamical properties of gases -----------

R_universal=8.31446261815324 # Universal gas constant, J.K-1.mol-1, should probably go elsewhere since it's a constant
    
def L_heat(switch): # Molar latent heat for gases considered, in J.mol-1
    if switch == 'H2O':
        return phys.water.L_vaporization*phys.water.MolecularWeight*1e-3 # Conversion from J.kg-1 to J.mol-1 (molecular weight is in g/mol in phys.f90)
    if switch == 'CH4':
        return phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
    if switch == 'CO2':
        return phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
    if switch == 'CO':
        return phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
    if switch == 'N2':
        return phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
    if switch == 'O2':
        return phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
    if switch == 'H2':
        return phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
    if switch == 'He':
        return phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
    if switch == 'NH3':
        return phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3

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
        return Tref/(1-(Tref*R_universal/L_heat('H2O'))*math.log(p/pref))
        #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
    if switch == 'CH4':
        Tref = 148.15 # K, arbitrary point (148.15K,esat(148.15K)=9.66bar) on the L/G coexistence curve of methane 
        pref = esat('CH4',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('CH4'))*math.log(p/pref))
    if switch == 'CO2':
        Tref = 253. # K, arbitrary point (253K,esat(253K)=20.9bar) on the coexistence curve of CO2 
        pref = esat('CO2',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('CO2'))*math.log(p/pref))
    if switch == 'CO':
        Tref = 100. # K, arbitrary point (100K,esat(100K)=4.6bar) on the coexistence curve of CO 
        pref = esat('CO',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('CO'))*math.log(p/pref))
    if switch == 'N2':
        Tref = 98.15 # K, arbitrary point (98.15K,esat(98.15K)=7.9bar) on the coexistence curve of N2 
        pref = esat('N2',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('N2'))*math.log(p/pref))
    if switch == 'O2':
        Tref = 123.15 # K, arbitrary point (123.15K,esat(123.15K)=21.9bar) on the coexistence curve of O2 
        pref = esat('O2',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('O2'))*math.log(p/pref))
    if switch == 'H2':
        Tref = 23.15 # K, arbitrary point (23.15K,esat(23.15K)=1.7bar) on the coexistence curve of H2 
        pref = esat('H2',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('H2'))*math.log(p/pref))
    if switch == 'He':
        Tref = 4.22 # K, boiling point of He at 1 atm 
        pref = 1e5 # esat('He',Tref) returns 45196 Pa, should return 1e5
        return Tref/(1-(Tref*R_universal/L_heat('He'))*math.log(p/pref))
    if switch == 'NH3':
        Tref = 273.15 # K, arbitrary point (273.15K,esat(273.15K)=8.6bar) on the coexistence curve of NH3 
        pref = esat('NH3',Tref)
        return Tref/(1-(Tref*R_universal/L_heat('NH3'))*math.log(p/pref))
    
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
    
def cpc(switch): # Molar heat capacities for condensed phases considered. J.K-1.mol-1. # ToDo: replace with liquid or solid phase values
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

def xv(switch,p,T): # Molar abundances of vapor phase/mole of gas mixture. p is the total pressure
    # if condensation occurs, the partial pressure pp=xv*p follows the Clausius-Clapeyron relation      
    if switch == 'H2O':       
        if T<Tdew('H2O',esat('H2O',T)):
            return esat('H2O',T)/p
        else:
            return 0.6
    if switch == 'CH4':
        if T<Tdew('CH4',esat('CH4',T)):          
            return esat('CH4',T)/p
        else:
            return 0.      
    if switch == 'CO2':
        if T<Tdew('CO2',esat('CO2',T)):           
            return esat('CO2',T)/p
        else:
            return 0.3
    if switch == 'CO':
        if T<Tdew('CO',esat('CO',T)):           
            return esat('CO',T)/p
        else:
            return 0.
    if switch == 'N2':
        if T<Tdew('N2',esat('N2',T)):        
            return esat('N2',T)/p
        else:
            return 0.6
    if switch == 'O2':
        if T<Tdew('O2',esat('O2',T)):          
            return esat('O2',T)/p
        else:
            return 0.
    if switch == 'H2':
        if T<Tdew('H2',esat('H2',T)):          
            return esat('H2',T)/p
        else:
            return 0.
    if switch == 'He':
        if T<Tdew('He',esat('He',T)):          
            return esat('He',T)/p
        else:
            return 0.
    if switch == 'NH3':
        if T<Tdew('NH3',esat('NH3',T)):          
            return esat('NH3',T)/p
        else:
            return 0.
    
def xc(switch): # Molar abundances of condensed phase/mole of gas mixture
    if switch == 'H2O':
        return 0.
    if switch == 'CH4':
        return 0.
    if switch == 'CO2':
        return 0.
    if switch == 'CO':
        return 0.
    if switch == 'N2':
        return 0.
    if switch == 'O2':
        return 0.
    if switch == 'H2':
        return 0.
    if switch == 'He':
        return 0.
    if switch == 'NH3':
        return 0.
    
    
#function that returns dT/dlnP
def General_moist_adiabat(T,lnP):
    # Add terms as they condense
    p = np.exp(lnP) # total pressure
    
    xd=1.-(xv('H2O',p,T)+\
           xv('CO2',p,T)+\
           xv('CH4',p,T)+\
           xv('CO',p,T) +\
           xv('N2',p,T) +\
           xv('O2',p,T) +\
           xv('H2',p,T) +\
           xv('He',p,T) +\
           xv('NH3',p,T)) # Molar abundance of non-condensible species/mole of gas mixture
    
    #xd=1.-(xv('H2O',p,T)+\
    #       xv('CO2',p,T)) # Molar abundance of non-condensible species/mole of gas mixture    
    
    cpd=1000.*phys.air.MolecularWeight*1.e-3 # molar heat capacity of the dry species
    #cpd=xv('N2',p,T)*cpv('N2') + xv('CO',p,T)*cpv('CO') + xv('CH4',p,T)*cpv('CH4') +\
    #                             xv('O2',p,T)*cpv('O2') + xv('H2',p,T)*cpv('H2') +\
    #                             xv('He',p,T)*cpv('He') + xv('NH3',p,T)*cpv('NH3') # molar heat capacity of the dry species
    
    # Dry version: xc=0, L terms=0
    if T>Tdew('H2O',p): 
        dTdlnP=R_universal*(xd+xv('H2O',p,T)+\
                               xv('CO2',p,T)+\
                               xv('CH4',p,T)+\
                               xv('CO',p,T) +\
                               xv('N2',p,T) +\
                               xv('O2',p,T) +\
                               xv('H2',p,T) +\
                               xv('He',p,T) +\
                               xv('NH3',p,T))/(xd*cpd+xv('H2O',p,T)*cpv('H2O')+\
                                                              xv('CO2',p,T)*cpv('CO2')+\
                                                              xv('CH4',p,T)*cpv('CH4')+\
                                                              xv('CO',p,T)*cpv('CO')+\
                                                              xv('N2',p,T)*cpv('N2')+\
                                                              xv('O2',p,T)*cpv('O2')+\
                                                              xv('H2',p,T)*cpv('H2')+\
                                                              xv('He',p,T)*cpv('He')+\
                                                              xv('NH3',p,T)*cpv('NH3'))
    
    if T<Tdew('H2O',p) and T>Tdew('CO2',p):
        xv_cpv=xv('H2O',p,T)*cpv('H2O') + xv('CO2',p,T)*cpv('CO2')+\
                                                xv('CH4',p,T)*cpv('CH4')+\
                                                xv('CO',p,T)*cpv('CO')+\
                                                xv('N2',p,T)*cpv('N2')+\
                                                xv('O2',p,T)*cpv('O2')+\
                                                xv('H2',p,T)*cpv('H2')+\
                                                xv('He',p,T)*cpv('He')+\
                                                xv('NH3',p,T)*cpv('NH3')
                                              
        first_term=(xv('H2O',p,T)/xd)*(L_heat('H2O')/(R_universal*T))**2
        quadr_term=((xv('H2O',p,T)/xd)*L_heat('H2O')/(R_universal*T))**2
        sum_abundances=xv('H2O',p,T) + xv('CO2',p,T)+\
                                             xv('CH4',p,T)+\
                                             xv('CO',p,T)+\
                                             xv('N2',p,T)+\
                                             xv('O2',p,T)+\
                                             xv('H2',p,T)+\
                                             xv('He',p,T)+\
                                             xv('NH3',p,T)
                                             
        sum_ratio_abundances=xv('H2O',p,T)/xd + xv('CO2',p,T)/xd+\
                                                      xv('CH4',p,T)/xd+\
                                                      xv('CO',p,T)/xd+\
                                                      xv('N2',p,T)/xd+\
                                                      xv('O2',p,T)/xd+\
                                                      xv('H2',p,T)/xd+\
                                                      xv('He',p,T)/xd+\
                                                      xv('NH3',p,T)/xd
        
        num=(1.+((xv('H2O',p,T)/xd)*L_heat('H2O')/(R_universal*T)))*T # Factor T to get dT/dlnP
        denom=(1./R_universal)*(xd*cpd+xv_cpv)/(xd+sum_abundances) + (first_term+quadr_term)/(1.+sum_ratio_abundances)
        dTdlnP=num/denom

    if T<Tdew('CO2',p): # and T>Tdew('N2'): to add once it works
        xv_cpv=xv('H2O',p,T)*cpv('H2O') + xv('CO2',p,T)*cpv('CO2')+\
                                                xv('CH4',p,T)*cpv('CH4')+\
                                                xv('CO',p,T)*cpv('CO')+\
                                                xv('N2',p,T)*cpv('N2')+\
                                                xv('O2',p,T)*cpv('O2')+\
                                                xv('H2',p,T)*cpv('H2')+\
                                                xv('He',p,T)*cpv('He')+\
                                                xv('NH3',p,T)*cpv('NH3')
                                              
        first_term=(xv('H2O',p,T)/xd)*(L_heat('H2O')/(R_universal*T))**2+\
                   (xv('CO2',p,T)/xd)*(L_heat('CO2')/(R_universal*T))**2
                   
        quadr_term=((xv('H2O',p,T)/xd)*L_heat('H2O')/(R_universal*T)+\
                    (xv('CO2',p,T)/xd)*L_heat('CO2')/(R_universal*T))**2
                    
        sum_abundances=xv('H2O',p,T) + xv('CO2',p,T)+\
                                             xv('CH4',p,T)+\
                                             xv('CO',p,T)+\
                                             xv('N2',p,T)+\
                                             xv('O2',p,T)+\
                                             xv('H2',p,T)+\
                                             xv('He',p,T)+\
                                             xv('NH3',p,T)
                                             
        sum_ratio_abundances=xv('H2O',p,T)/xd + xv('CO2',p,T)/xd+\
                                                      xv('CH4',p,T)/xd+\
                                                      xv('CO',p,T)/xd+\
                                                      xv('N2',p,T)/xd+\
                                                      xv('O2',p,T)/xd+\
                                                      xv('H2',p,T)/xd+\
                                                      xv('He',p,T)/xd+\
                                                      xv('NH3',p,T)/xd
                                            
        num=(1.+((xv('H2O',p,T)/xd)*L_heat('H2O')/(R_universal*T) + (xv('CO2',p,T)/xd)*L_heat('CO2')/(R_universal*T)))*T
        denom=(1./R_universal)*(xd*cpd+xv_cpv)/(xd+sum_abundances) + (first_term+quadr_term)/(1.+sum_ratio_abundances)
        dTdlnP=num/denom

    # ToDo: add other versions with more condensibles
       
    return dTdlnP
       
#lnP points
nlev = 100         # number of pressure levels
ps = 1.e5       # Surface pressure in Pa
ptop = 1e-5*ps   # Pressure at the TOA in Pa
p = np.exp(np.linspace(math.log(ps),math.log(ptop),nlev)) # Logarithmic pressure grid scale

# Ground temperature
Tg=300.

#solve ODE
General_moist_adiabat=odeint(General_moist_adiabat,Tg,np.log(p)) 

#plot results
H2O_saturation=[Tdew('H2O',pp) for pp in p] # just to plot it

plt.semilogy(General_moist_adiabat,p)
#plt.semilogy(H2O_saturation,p, label='H2O saturation vapor pressure curve')
plt.gca().invert_yaxis()
plt.xlabel('T(lnP)')
plt.ylabel('P')
plt.legend()
plt.show()