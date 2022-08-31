#!/usr/bin/env python
# coding: utf-8

# In[49]:


import phys 
import ClimateUtilities
import numpy as np
import math
import matplotlib.pyplot as plt
#%%
# Code for integrating simple case of no-condensate-retention adiabat

'''
Hi Ray. This code produces the adiabat I sent you. I think all of the import
calls are things you should have easy access to. 

With the current settings, it should just run and produce the plot of the adiabat.

Sorry for the annoying "cpv", "L_heat", and "p_sat" functions, had to yank those 
from other pieces of code to avoid sending you a bundle of different scripts.

'''
#get_ipython().run_line_magic('matplotlib', 'notebook')


# In[50]:


R_gas = 8.31446261815324

#Vapor specific heat function

def cpv( vol, tmp, cp_mode = "constant" ):

    
    # cp_mode = "T-dependent" # NIST Chemistry WebBook

    if cp_mode == "T-dependent":

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2
        if vol == "H2O":
            # Temperature (K) 500. - 1700.    1700. - 6000.
            A = [ 30.09200  ,  41.96426  ]
            B = [ 6.832514  ,  8.622053  ]
            C = [ 6.793435  ,  -1.499780 ] 
            D = [ -2.534480 ,  0.098119  ]
            E = [ 0.082139  ,  -11.15764 ] 
            F = [ -250.8810 ,  -272.1797 ] 
            G = [ 223.3967  ,  219.7809  ]
            H = [ -241.8264 ,  -241.8264 ] 
            if tmp <= 1700:
                cp_idx = 0
            if tmp > 1700:
                cp_idx = 1
            tmp = np.max([tmp, 500]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/CO2/c2-1-3
        if vol == "CO2":
            # Temperature (K) 298. - 1200.    1200. - 6000.
            A = [ 24.99735  , 58.16639  ]
            B = [ 55.18696  , 2.720074  ]
            C = [ -33.69137 , -0.492289 ] 
            D = [ 7.948387  , 0.038844  ]
            E = [ -0.136638 , -6.447293 ] 
            F = [ -403.6075 , -425.9186 ] 
            G = [ 228.2431  , 263.6125  ]
            H = [ -393.5224 , -393.5224 ] 
            if tmp <= 1200:
                cp_idx = 0
            if tmp > 1200:
                cp_idx = 1
            tmp = np.max([tmp, 298]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2/h1H
        if vol == "H2":
            # Temperature (K) 298. - 1000.    1000. - 2500.   2500. - 6000.
            A = [ 33.066178  , 18.563083  , 43.413560  ]      
            B = [ -11.363417 , 12.257357  , -4.293079  ]      
            C = [ 11.432816  , -2.859786  , 1.272428   ]     
            D = [ -2.772874  , 0.268238   , -0.096876  ]      
            E = [ -0.158558  , 1.977990   , -20.533862 ]       
            F = [ -9.980797  , -1.147438  , -38.515158 ]       
            G = [ 172.707974 , 156.288133 , 162.081354 ]       
            H = [ 0.0        , 0.0        , 0.0        ] 
            if tmp <= 1000:
                cp_idx = 0
            if tmp > 1000 and tmp <= 2500:
                cp_idx = 1
            if tmp > 2500:
                cp_idx = 2
            tmp = np.max([tmp, 298]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/N2/c1-2
        if vol == "N2":
            # Temperature (K) 100. - 500. 500. - 2000.    2000. - 6000.
            A = [ 28.98641  ,  19.50583  ,  35.51872  ]    
            B = [ 1.853978  ,  19.88705  ,  1.128728  ]    
            C = [ -9.647459 ,  -8.598535 ,  -0.196103 ]     
            D = [ 16.63537  ,  1.369784  ,  0.014662  ]    
            E = [ 0.000117  ,  0.527601  ,  -4.553760 ]     
            F = [ -8.671914 ,  -4.935202 ,  -18.97091 ]     
            G = [ 226.4168  ,  212.3900  ,  224.9810  ]    
            H = [ 0.0       ,  0.0       ,  0.0       ]
            if tmp <= 500:
                cp_idx = 0
            if tmp > 500 and tmp <= 2000:
                cp_idx = 1
            if tmp > 2000:
                cp_idx = 2
            tmp = np.max([tmp, 100]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/CH4/h1H4
        if vol == "CH4":
            # Temperature (K) 298. - 1300.    1300. - 6000.
            A = [ -0.703029 ,  85.81217  ]  
            B = [ 108.4773  ,  11.26467  ]  
            C = [ -42.52157 ,  -2.114146 ]   
            D = [ 5.862788  ,  0.138190  ]  
            E = [ 0.678565  ,  -26.42221 ]   
            F = [ -76.84376 ,  -153.5327 ]   
            G = [ 158.7163  ,  224.4143  ]  
            H = [ -74.87310 ,  -74.87310 ]    
            if tmp <= 1300:
                cp_idx = 0
            if tmp > 1300:
                cp_idx = 1
            tmp = np.max([tmp, 298]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/CO/c1-2
        if vol == "CO":
            # Temperature (K) 298. - 1300.    1300. - 6000.
            A = [ 25.56759  ,  35.15070  ]   
            B = [ 6.096130  ,  1.300095  ]  
            C = [ 4.054656  ,  -0.205921 ]   
            D = [ -2.671301 ,  0.013550  ]  
            E = [ 0.131021  ,  -3.282780 ]   
            F = [ -118.0089 ,  -127.8375 ]   
            G = [ 227.3665  ,  231.7120  ]  
            H = [ -110.5271 ,  -110.5271 ]    
            if tmp <= 1300:
                cp_idx = 0
            if tmp > 1300:
                cp_idx = 1
            tmp = np.max([tmp, 298]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/O2/c1-2
        if vol == "O2":
            # Temperature (K) 100. - 700. 700. - 2000.    2000. - 6000.
            A = [ 31.32234  ,  30.03235  ,  20.91111  ]       
            B = [ -20.23531 ,  8.772972  ,  10.72071  ]       
            C = [ 57.86644  ,  -3.988133 ,  -2.020498 ]        
            D = [ -36.50624 ,  0.788313  ,  0.146449  ]       
            E = [ -0.007374 ,  -0.741599 ,  9.245722  ]       
            F = [ -8.903471 ,  -11.32468 ,  5.337651  ]       
            G = [ 246.7945  ,  236.1663  ,  237.6185  ]       
            H = [ 0.0       ,  0.0       ,  0.0       ]  
            if tmp <= 700:
                cp_idx = 0
            if tmp > 700 and tmp <= 2000:
                cp_idx = 1
            if tmp > 2000:
                cp_idx = 2
            tmp = np.max([tmp, 100]) # Fit validity
        
        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/He
        if vol == "He":
            # Temperature (K) 298. - 6000.
            A = [ 20.78603     ]
            B = [ 4.850638e-10 ]    
            C = [ -1.582916e-10]     
            D = [ 1.525102e-11 ]    
            E = [ 3.196347e-11 ]    
            F = [ -6.197341    ] 
            G = [ 151.3064     ]
            H = [ 0.000000     ] 
            cp_idx = 0
            tmp = np.max([tmp, 298]) # Fit validity

        # https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H3N/h1H3
        if vol == "NH3":
            # Temperature (K) 298. - 1400.    1400. - 6000.
            A = [ 19.99563  ,  52.02427  ]   
            B = [ 49.77119  ,  18.48801  ]   
            C = [ -15.37599 ,  -3.765128 ]    
            D = [ 1.921168  ,  0.248541  ]   
            E = [ 0.189174  ,  -12.45799 ]    
            F = [ -53.30667 ,  -85.53895 ]    
            G = [ 203.8591  ,  223.8022  ]   
            H = [ -45.89806 ,  -45.89806 ]    
            if tmp <= 1400:
                cp_idx = 0
            if tmp > 1400:
                cp_idx = 1
            tmp = np.max([tmp, 298]) # Fit validity
            
        t = tmp/1000.
        cp = A[cp_idx] + B[cp_idx]*t + C[cp_idx]*t**2. + D[cp_idx]*t**3. + E[cp_idx]/t**2.

        return cp # J mol-1 K-1
    
    if cp_mode == "constant":
        if vol == 'H2O':
            cp = phys.water.cp*phys.water.MolecularWeight*1e-3
        if vol == 'CH4':
            cp = phys.methane.cp*phys.methane.MolecularWeight*1e-3
        if vol == 'CO2':
            cp = phys.co2.cp*phys.co2.MolecularWeight*1e-3
        if vol == 'CO':
            cp = phys.co.cp*phys.co.MolecularWeight*1e-3
        if vol == 'N2':
            cp = phys.N2.cp*phys.N2.MolecularWeight*1e-3
        if vol == 'O2':
            cp = phys.o2.cp*phys.o2.MolecularWeight*1e-3
        if vol == 'H2':
            cp = phys.h2.cp*phys.h2.MolecularWeight*1e-3
        if vol == 'He':
            cp = phys.he.cp*phys.he.MolecularWeight*1e-3
        if vol == 'NH3':
            cp = phys.nh3.cp*phys.nh3.MolecularWeight*1e-3   

        return cp # J mol-1 K-1 
#saturation vapor pressure function
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
        e = phys.satvps_function(phys.N2)
    if switch == 'O2':
        e = phys.satvps_function(phys.o2)
    if switch == 'H2':
        e = phys.satvps_function(phys.h2)
    if switch == 'He':
        e = phys.satvps_function(phys.he)
    if switch == 'NH3':
        e = phys.satvps_function(phys.nh3)
    
    # Return saturation vapor pressure
    #return float(f'{e(T):.2f}') #RTP: Why do the return this way? 
    return e(T)
#Latent Heat function
def L_heat(switch, T):

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

    # # No latent heat contribution in moist adiabat if below p_sat
    # if P < p_sat(switch, T):
    # # if T > Tdew(switch, psat):
    #     L_heat = 0.

    return L_heat  


# In[62]:


#Simplified Moist Slope function. No condensate retention. Takes lnP, lnT, and
# a list of volatile fractions in molar units (vol_list)
def moist_slope_no_atm_no_cond(lnP, lnT, vol_list):
    # T instead lnT
    tmp = math.exp(lnT)

    # Sum terms in equation
    num_sum     = 0.
    denom_sum1  = 0. 
    denom_sum2  = 0. 
    denom_sum3  = 0.
    cp          = 0.
    xd          = 0.
    xv          = 0.
    #Check if volatile is subsaturated, or (super)saturated based on p_sat
    for vol in vol_list.keys():
        p_vol=np.exp(lnP)*vol_list[vol]       
        if np.isclose(p_vol, p_sat(vol,tmp)):
            xv += vol_list[vol]
            #print(vol + 'is saturated')
        elif p_vol < p_sat(vol, tmp):            
            xd += vol_list[vol]
            #print(vol + 'is subsaturated')
        elif p_vol > p_sat(vol,tmp):
            xv += vol_list[vol]
            #print(vol + 'is supersaturated')
            #print('Warning: volatile ' + vol + ' is supersaturated. Psat=%.3f'%p_sat(vol,tmp)+', Pvol=%.3f'%p_vol)
    
    # Calculate sums over volatiles
    for vol in vol_list.keys(): 
        p_vol=np.exp(lnP)*vol_list[vol]
        # Coefficients
        eta_vol     = vol_list[vol] / xd
        
        if np.isclose(p_vol,p_sat(vol,tmp)) or p_vol > p_sat(vol,tmp) or p_vol == p_sat(vol,tmp):
            beta_vol    = L_heat(vol, tmp) / (R_gas * tmp) #RTP

        # Beta terms zero if below saturation vapor pressure
        else: 
            beta_vol = 0.
        # Sum in numerator
        num_sum     += eta_vol * beta_vol
        #print(num_sum)
        #print(beta_vol)
        # Sums in denominator
        denom_sum1  += eta_vol * (beta_vol**2.)
        denom_sum3  += eta_vol
        cp    += vol_list[vol] * cpv(vol, tmp)
        
    cp = cp / ( xd + xv )
    
    # Sum 2 in denominator  
    denom_sum2  = num_sum ** 2.

    # Collect terms
    numerator   = 1. + num_sum
    denominator = (cp / R_gas) + (denom_sum1 + denom_sum2) / (1. + denom_sum3) #RTP R_gas

    # dlnT/dlnP
    dlnTdlnP = numerator / denominator

    # Moist adiabat slope
    return dlnTdlnP


# In[69]:


def dlnT_dlnP_d(lnP, lnT, vol_list):
    # T instead lnT
    tmp = math.exp(lnT)

    # Sum terms in equation
    num_sum     = 0.
    denom_sum  = 0. 
    cp          = 0.
    xd          = 0.
    xv          = 0.
    #Check if volatile is subsaturated, or (super)saturated based on p_sat
    for vol in vol_list.keys():
        p_vol=np.exp(lnP)*vol_list[vol]       
        if np.isclose(p_vol, p_sat(vol,tmp)):
            xv += vol_list[vol]
        elif p_vol < p_sat(vol, tmp):            
            xd += vol_list[vol]
        elif p_vol > p_sat(vol,tmp):
            xv += vol_list[vol]
    
    # Calculate sums over volatiles
    for vol in vol_list.keys(): 
        p_vol=np.exp(lnP)*vol_list[vol]
        # Coefficients
        eta_vol     = vol_list[vol] / xd
        #print(eta_vol)
        # sums for saturated comps
        if np.isclose(p_vol,p_sat(vol,tmp)) or p_vol > p_sat(vol,tmp):
            L = L_heat(vol,tmp)
            #beta_vol    = L_heat(vol, tmp) / (R_gas * tmp) #RTP
            # Sum in numerator
            num_sum     += eta_vol * L / tmp
            # Sum in denominator
            denom_sum  += eta_vol * (cpv(vol, tmp) - L/tmp + L**2/(R_gas*tmp**2))
            
        # sums for subsaturated comps
        elif p_vol < p_sat(vol, tmp): 
            #print(eta_vol)
            denom_sum += cpv(vol,tmp) * eta_vol #The eta_vol (x_vol/x_d) is there so that if there are multiple dry components the average dry specific heat is used
        
    
        
    
    # Collect terms
    numerator   = R_gas + num_sum
    denominator = denom_sum

    # dlnT/dlnPd
    dlnTdlnPd = numerator / denominator

    # Moist adiabat slope
    return dlnTdlnPd


# In[70]:


def dlnT_dlnP_from_product(lnP, lnT, vol_list):
    # T instead lnT
    tmp = math.exp(lnT)

    # Sum terms in equation
    num_sum     = 0.
    denom_sum  = 0. 
    cp          = 0.
    xd          = 0.
    xv          = 0.
    #Check if volatile is subsaturated, or (super)saturated based on p_sat
    for vol in vol_list.keys():
        p_vol=np.exp(lnP)*vol_list[vol]       
        if np.isclose(p_vol, p_sat(vol,tmp)):
            xv += vol_list[vol]
            #print(vol + 'is saturated')
        elif p_vol < p_sat(vol, tmp):            
            xd += vol_list[vol]
            #print(vol + 'is subsaturated')
        elif p_vol > p_sat(vol,tmp):
            xv += vol_list[vol]
            #print(vol + 'is supersaturated')
            #print('Warning: volatile ' + vol + ' is supersaturated. Psat=%.3f'%p_sat(vol,tmp)+', Pvol=%.3f'%p_vol)
    
    # Calculate sums over volatiles
    for vol in vol_list.keys(): 
        p_vol=np.exp(lnP)*vol_list[vol]
        # Coefficients
        eta_vol     = vol_list[vol] / xd
        #print(eta_vol)
        # sums for saturated comps
        if np.isclose(p_vol,p_sat(vol,tmp)) or p_vol > p_sat(vol,tmp):
            L = L_heat(vol,tmp)
            
            # Sum in numerator
            num_sum     += eta_vol
            # Sum in denominator
            denom_sum  += L/R_gas/tmp * eta_vol * dlnT_dlnP_d(lnP,lnT,vol_list)
            
        
        
    
        
    
    # Collect terms
    numerator   = 1 + num_sum
    denominator = 1 + denom_sum
    dlnPd_dlnP = numerator / denominator
    # dlnT/dlnPd
    dlnTdlnP = dlnPd_dlnP * dlnT_dlnP_d(lnP,lnT,vol_list)

    # Moist adiabat slope
    return dlnTdlnP


# ### Code for RTP adiabat slope function
# 
# To convert from ClimateBook textbook formula for $d\ln T/d\ln p_a$ to $d\ln T/d\ln p$, divide by
# 
# $$ \frac{d\ln p}{d\ln p_a} = \frac{p_a}{p} (1 + \frac{L}{R_c T}\frac{p_c}{p_a} \frac{d\ln T}{d\ln p_a} )$$
# 
# which is derived using $p = p_a + p_{sat}(T)$ and Clausius-Clapeyron.
# 
# Suggestion to RJ:  Put in my formula for $d\ln T/d\ln p_a$, and carry out the algebra to get the expression for $d\ln T/d\ln p$, then compare with the RJ adiabat formula with zero retained condensate.  
# 
# Side note, after conversation with David Romps: Does the RJ adiabat formula equal the Li/Ingersoll slope formula when condensate is set to zero?  David thinks that the derivation of the pseudo adiabat should give a different slope than setting retained condensate to zero in Li/Ingersoll, unless Li/Ingersoll made some special additional assumptions.  Worth checking out after we sort out the basic discrepancy. 

# In[71]:


#RTP functions
condensible = phys.H2O
noncondensible = phys.N2
satvp = phys.satvps_function(condensible,'liquid')#phys.satvps_function(condensible) -- Changed for 2d ed.
#Set up thermodynamic constants
eps = condensible.MolecularWeight/noncondensible.MolecularWeight
L = condensible.L_vaporization
Ra = noncondensible.R
Rc = condensible.R
cpa = noncondensible.cp
cpc = condensible.cp
def slopeRTP(logpaIn,logTIn):
    pa1 = math.exp(logpaIn)
    T1 = math.exp(logTIn)
    rsat1 = eps*(satvp(T1)/pa1) #mixing ratio
    num = (1. + (L/(Ra*T1))*rsat1)*Ra
    den = cpa + (cpc + (L/(Rc*T1) - 1.)*(L/T1))*rsat1
    return num/den
def slopeRTP_ptot(logp,logT):
    p,T = math.exp(logp),math.exp(logT)
    pa = p - satvp(T)
    logpa = math.log(pa)
    b = L/(Rc*T)
    slope0 = slopeRTP(logpa,logT)
    slope = slope0*(p/pa)/(1.+ b*(satvp(T)/pa)*slope0)
    return slope


# In[72]:


#pN2,T = 1.e5,350.
#slopeRTP(np.log(pN2),np.log(T)), slopeRTP_ptot(np.log(pN2+satvp(T)),np.log(T)),(phys.water.R*T)/phys.water.L_vaporization


# In[78]:

plt.figure()
pN2 = 1.e5
vol_list = {}
temp_arr = np.arange(274, 500)
p_arr = pN2+np.vectorize(satvp)(temp_arr)
RJ_slope_arr = []
RJ_old_slope_arr = []
Ray_slope_arr = []
for i in range(len(temp_arr)):
    press = p_arr[i]
    Temp = temp_arr[i]
    vol_list['N2'] = pN2/press
    vol_list['H2O'] = satvp(Temp)/press
    RJ_slope_arr.append(dlnT_dlnP_from_product(np.log(press),np.log(Temp),vol_list))
    RJ_old_slope_arr.append(moist_slope_no_atm_no_cond(np.log(press), np.log(Temp), vol_list))
    Ray_slope_arr.append(slopeRTP_ptot(np.log(press),np.log(Temp)))
    
#moist_slope_no_atm_no_cond(np.log(p), np.log(T), vol_list)
#print(dlnT_dlnP_from_product(np.log(p),np.log(T),vol_list),slopeRTP_ptot(np.log(pN2+satvp(T)),np.log(T)))


plt.plot(temp_arr,RJ_slope_arr,label='Slope from product')
plt.plot(temp_arr,RJ_old_slope_arr,label='Slope in paper\'s form')
plt.plot(temp_arr,Ray_slope_arr,label='Ray')
plt.legend(loc='best')
plt.ylabel(r'$\frac{dlnT}{dlnP}$')
plt.xlabel('Temperature [K]')
plt.title('Adiabat slopes vs. temperature for atmosphere with 1 bar \n N2 background and saturated H2O')
# ## Integrate, and make plots of RJ adiabat

# In[76]:


#%%
# Surface conditions:
pN2_surf = 1e5
T_surf                  = 275.         # K
P_surf                  = p_sat('H2O',T_surf) + pN2_surf      # Pa
#integration step
step = -0.01
# Top of atmosphere pressure
ptop = 100 #Pa
Tmin = 75 #K
# Volatile molar concentrations: ! should sum to one !
vol_list = { 
              "H2O" : p_sat('H2O',T_surf)/P_surf,        
              "N2"  : pN2_surf/P_surf,      
              }
#lists for populating during integration
moist_tuple = []
pressure_list = []
temp_list = []
pressure_list.append(P_surf)
temp_list.append(T_surf)

#Create integrator instance
int_slope = ClimateUtilities.integrator(dlnT_dlnP_from_product,np.log(P_surf), np.log(T_surf), step)

# Feed it the volatile fractions
int_slope.setParams(vol_list)


while pressure_list[-1] > ptop and temp_list[-1]>Tmin:
    #Integrate!
    moist_tuple.append(int_slope.next())
    #Feed the new values to the respective lists
    pressure_list.append(np.exp(int_slope.x))
    temp_list.append(np.exp(int_slope.y))
    #New p_h2o !ASSUMES WATER IS SATURATED FOR THE PURPOSES OF THIS SIMPLE BIT OF CODE
    p_h2o = p_sat('H2O',temp_list[-1])
    #New p_n2 ! ASSUMES N2 SUBSATURATED FOR THIS SIMPLE CODE
    p_n2 = pressure_list[-1] - p_h2o
    # New vol_list with new mole fractions
    vol_list = { 
              "H2O" : p_h2o/pressure_list[-1],        
              "N2"  : p_n2/pressure_list[-1],       
              }
    # set parameters in int_slope to new vol_list (different results if I don't do this, doesn't seem to autoupdate)
    int_slope.setParams(vol_list)
    #print(pressure_list[-1],vol_list['N2'],vol_list['H2O'])

#%%
plt.figure()
#%%
plt.plot(temp_list,pressure_list, label = 'New adiabat formulation')
#%%
plt.yscale('log')
plt.gca().invert_yaxis()
plt.title('Temperature profile for 1 bar N2 + saturated H2O')
plt.ylabel('Pressure [Pa]')
plt.xlabel('Temperature [K]')


# In[ ]:
P_surf                  = p_sat('H2O',T_surf) + pN2_surf      # Pa
#integration step
step = -0.01
# Top of atmosphere pressure
ptop = 100 #Pa
# Volatile molar concentrations: ! should sum to one !
vol_list = { 
              "H2O" : p_sat('H2O',T_surf)/P_surf,        
              "N2"  : pN2_surf/P_surf,      
              }
#lists for populating during integration
moist_tuple = []
pressure_list = []
temp_list = []
pressure_list.append(P_surf)
temp_list.append(T_surf)
#Create integrator instance
int_slope = ClimateUtilities.integrator(moist_slope_no_atm_no_cond,np.log(P_surf), np.log(T_surf), step)

# Feed it the volatile fractions
int_slope.setParams(vol_list)


while pressure_list[-1] > ptop and temp_list[-1]>Tmin:
    #Integrate!
    moist_tuple.append(int_slope.next())
    #Feed the new values to the respective lists
    pressure_list.append(np.exp(int_slope.x))
    temp_list.append(np.exp(int_slope.y))
    #New p_h2o !ASSUMES WATER IS SATURATED FOR THE PURPOSES OF THIS SIMPLE BIT OF CODE
    p_h2o = p_sat('H2O',temp_list[-1])
    #New p_n2 ! ASSUMES N2 SUBSATURATED FOR THIS SIMPLE CODE
    p_n2 = pressure_list[-1] - p_h2o
    # New vol_list with new mole fractions
    vol_list = { 
              "H2O" : p_h2o/pressure_list[-1],        
              "N2"  : p_n2/pressure_list[-1],       
              }
    # set parameters in int_slope to new vol_list (different results if I don't do this, doesn't seem to autoupdate)
    int_slope.setParams(vol_list)
    #print(pressure_list[-1],vol_list['N2'],vol_list['H2O'])


plt.plot(temp_list,pressure_list, label = 'Paper\'s adiabat formulation')
#%%
P_surf                  = p_sat('H2O',T_surf) + pN2_surf      # Pa
#integration step
step = -0.01
# Top of atmosphere pressure
ptop = 100 #Pa
# Volatile molar concentrations: ! should sum to one !
vol_list = { 
              "H2O" : p_sat('H2O',T_surf)/P_surf,        
              "N2"  : pN2_surf/P_surf,      
              }
#lists for populating during integration
moist_tuple = []
pressure_list = []
temp_list = []
pressure_list.append(P_surf)
temp_list.append(T_surf)

#Create integrator instance
int_slope = ClimateUtilities.integrator(slopeRTP_ptot,np.log(P_surf), np.log(T_surf), step)

# Feed it the volatile fractions
int_slope.setParams(vol_list)


while pressure_list[-1] > ptop and temp_list[-1]>Tmin:
    #Integrate!
    moist_tuple.append(int_slope.next())
    #Feed the new values to the respective lists
    pressure_list.append(np.exp(int_slope.x))
    temp_list.append(np.exp(int_slope.y))
    #New p_h2o !ASSUMES WATER IS SATURATED FOR THE PURPOSES OF THIS SIMPLE BIT OF CODE
    p_h2o = p_sat('H2O',temp_list[-1])
    #New p_n2 ! ASSUMES N2 SUBSATURATED FOR THIS SIMPLE CODE
    p_n2 = pressure_list[-1] - p_h2o
    # New vol_list with new mole fractions
    vol_list = { 
              "H2O" : p_h2o/pressure_list[-1],        
              "N2"  : p_n2/pressure_list[-1],       
              }
    # set parameters in int_slope to new vol_list (different results if I don't do this, doesn't seem to autoupdate)
    int_slope.setParams(vol_list)
    #print(pressure_list[-1],vol_list['N2'],vol_list['H2O'])


plt.plot(temp_list,pressure_list, label = 'Ray adiabat')
plt.legend(loc='best')

# In[ ]:





# In[ ]:




