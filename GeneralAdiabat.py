'''
Created 01/12/19

@authors: 
Ryan Boukrouche (RB)
Tim Lichtenberg (TL)
RJ Graham (RJ)
This file builds a self-consistent temperature profile from either Li, Ingersoll & Oyafuso 2018 or 
Graham, Boukrouche, Lichtenberg, Pierrehumbert 2020, from the ground up using the Runge-Kutta 4 scheme 
from ClimateUtilities.py.
'''

import time
import numpy as np
import math
# from scipy.integrate import odeint
# from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# import seaborn as sns
import copy
from matplotlib import cm

# Coupler-specific modules
try:
    from ClimateUtilities import *
    from atmosphere_column import atmos
    import phys
except:
    from atm_rad_conv.ClimateUtilities import *
    from atm_rad_conv.atmosphere_column import atmos
    import atm_rad_conv.phys as phys

# Color definitions: 
# https://www.codecademy.com/articles/seaborn-design-ii
# https://python-graph-gallery.com/python-colors/
# https://matplotlib.org/tutorials/colors/colormaps.html
# https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
no_colors   = 7
vol_colors = {
    "H2O"            : cm.get_cmap('PuBu', no_colors)(range(no_colors)),
    "CO2"            : cm.get_cmap("Reds", no_colors)(range(no_colors)),
    "H2"             : cm.get_cmap("Greens", no_colors)(range(no_colors)),
    "N2"             : cm.get_cmap("Purples", no_colors)(range(no_colors)),
    "O2"             : cm.get_cmap("Wistia", no_colors+2)(range(no_colors+2)),
    "CH4"            : cm.get_cmap("RdPu", no_colors)(range(no_colors)),
    "CO"             : cm.get_cmap("pink_r", no_colors)(range(no_colors)),
    "S"              : cm.get_cmap("YlOrBr", no_colors)(range(no_colors)),
    "He"             : cm.get_cmap("Greys", no_colors)(range(no_colors)),
    "NH3"            : cm.get_cmap("cool", no_colors)(range(no_colors)),
    "mixtures"       : cm.get_cmap("Set3", 9)(range(no_colors)),
    "H2O-CO2"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "CO2-H2O"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "H2O-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-H2O"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-CO"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "CO-H2"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "H2-CO2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "CO2-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "H2-CH4"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "CH4-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "H2-N2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "N2-H2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "CO2-N2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    "N2-CO2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    # "H2O"            : sns.color_palette("PuBu", no_colors),
    # "CO2"            : sns.color_palette("Reds", no_colors),
    # "H2"             : sns.color_palette("Greens", no_colors),
    # "N2"             : sns.color_palette("Purples", no_colors), # sns.cubehelix_palette(7)
    # "O2"             : sns.light_palette("darkturquoise", no_colors),
    # "CH4"            : sns.color_palette("RdPu", no_colors),
    # "CO"             : sns.light_palette("#731d1d", no_colors),
    # "S"              : sns.light_palette("#EBB434", no_colors),
    # "He"             : sns.color_palette("Greys", no_colors),
    # "NH3"            : sns.light_palette("teal", no_colors),
    # "mixtures"       : sns.color_palette("Set2", 11),
    # "H2O-CO2"        : sns.color_palette("Set2", 11)[9],
    # "CO2-H2O"        : sns.color_palette("Set2", 11)[9],
    # "H2O-H2"         : sns.color_palette("Set2", 11)[3],
    # "H2-H2O"         : sns.color_palette("Set2", 11)[3],
    # "H2-CO"          : sns.color_palette("Set2", 11)[7],
    # "CO-H2"          : sns.color_palette("Set2", 11)[7],
    # "H2-CO2"         : sns.color_palette("Set2", 11)[8],
    # "CO2-H2"         : sns.color_palette("Set2", 11)[8],
    # "H2-CH4"         : sns.color_palette("Set2", 11)[2],
    # "CH4-H2"         : sns.color_palette("Set2", 11)[2],
    # "H2-N2"          : sns.color_palette("Set2", 11)[4],
    # "N2-H2"          : sns.color_palette("Set2", 11)[4],
    # "CO2-N2"         : sns.color_palette("Set2", 11)[5],
    # "N2-CO2"         : sns.color_palette("Set2", 11)[5],
    "black_1"        : "#000000",
    "black_2"        : "#323232",
    "black_3"        : "#7f7f7f",
    "qgray"          : "#768E95",
    "qgray2"         : "#888888",
    "qblue"          : "#4283A9", # http://www.color-hex.com/color/4283a9
    "qgreen"         : "#62B4A9", # http://www.color-hex.com/color/62b4a9
    "qred"           : "#E6767A",
    "qturq"          : "#2EC0D1",
    "qorange"        : "#ff7f0e",
    "qmagenta"       : "#9A607F",
    "qyellow"        : "#EBB434",
    "qgray_dark"     : "#465559",
    "qblue_dark"     : "#274e65",
    "qgreen_dark"    : "#3a6c65",
    "qred_dark"      : "#b85e61",
    "qturq_dark"     : "#2499a7",
    "qmagenta_dark"  : "#4d303f",
    "qyellow_dark"   : "#a47d24",
    "qgray_light"    : "#acbbbf",
    "qblue_light"    : "#8db4cb",
    "qgreen_light"   : "#a0d2cb",
    "qred_light"     : "#eb9194",
    "qturq_light"    : "#57ccda",
    "qmagenta_light" : "#c29fb2",
    "qyellow_light"  : "#f1ca70",
}

# Volatile Latex names
vol_latex = {
    "H2O"     : r"H$_2$O",
    "CO2"     : r"CO$_2$",
    "H2"      : r"H$_2$" ,
    "CH4"     : r"CH$_4$",
    "CO"      : r"CO",
    "N2"      : r"N$_2$",
    "S"       : r"S",
    "O2"      : r"O$_2$",
    "He"      : r"He",
    "NH3"     : r"NH$_3$",
    "H2O-CO2" : r"H$_2$O–CO$_2$",
    "H2O-H2"  : r"H$_2$O–H$_2$",
    "H2O-CO"  : r"H$_2$O–CO",
    "H2O-CH4" : r"H$_2$O–CH$_4$",
    "H2O-N2"  : r"H$_2$O–N$_2$",
    "H2O-O2"  : r"H$_2$O–O$_2$",
    "H2-H2O"  : r"H$_2$–H$_2$O",
    "H2-CO"   : r"H$_2$–CO",
    "H2-CH4"  : r"H$_2$–CH$_4$",
    "H2-CO2"  : r"H$_2$–CO$_2$",
    "H2-N2"   : r"H$_2$–N$_2$",
    "H2-O2"   : r"H$_2$-O$_2$",
    "CO2-N2"  : r"CO$_2$–N$_2$",
    "CO2-H2O" : r"CO$_2$–H$_2$O",
    "CO2-CO"  : r"CO$_2$–CO",
    "CO2-CH4"  : r"CO$_2$–CH$_4$",
    "CO2-O2"  : r"CO$_2$–O$_2$",
    "CO2-H2"  : r"CO$_2$–H$_2$",
    "CO-H2O" : r"CO–H$_2$O",
    "CO-CO2" : r"CO–CO$_2$",
    "CO-H2"  : r"CO–H$_2$",
    "CO-CH4" : r"CO–CH$_4$",
    "CO-N2"  : r"CO–N$_2$",
    "CO-O2"  : r"CO–O$_2$",
    "CH4-H2O" : r"CH$_4$–H$_2$O",
    "CH4-CO2" : r"CH$_4$–CO$_2$",
    "CH4-H2"  : r"CH$_4$–H$_2$",
    "CH4-CO"  : r"CH$_4$–CO",
    "CH4-CH4" : r"CH$_4$–CH$_4$",
    "CH4-N2"  : r"CH$_4$–N$_2$",
    "CH4-O2"  : r"CH$_4$–O$_2$",
    "N2-H2O" : r"N$_2$–H$_2$O",
    "N2-CO2" : r"N$_2$–CO$_2$",
    "N2-H2"  : r"N$_2$–H$_2$",
    "N2-CO"  : r"N$_2$–CO",
    "N2-CH4" : r"N$_2$–CH$_4$",
    "N2-N2"  : r"N$_2$–N$_2$",
    "N2-O2"  : r"N$_2$–O$_2$",
    "O2-H2O" : r"O$_2$–H$_2$O",
    "O2-CO2" : r"O$_2$–CO$_2$",
    "O2-H2"  : r"O$_2$–H$_2$",
    "O2-CO"  : r"O$_2$–CO",
    "O2-CH4" : r"O$_2$–CH$_4$",
    "O2-N2"  : r"O$_2$–N$_2$",
    "O2-O2"  : r"O$_2$–O$_2$",
}

molar_mass      = {
          "H2O" : 0.01801528,           # kg mol−1
          "CO2" : 0.04401,              # kg mol−1
          "H2"  : 0.00201588,           # kg mol−1
          "CH4" : 0.01604,              # kg mol−1
          "CO"  : 0.02801,              # kg mol−1
          "N2"  : 0.028014,             # kg mol−1
          "O2"  : 0.031999,             # kg mol−1
          "SO2" : 0.064066,             # kg mol−1
          "H2S" : 0.0341,               # kg mol−1 
          "H"   : 0.001008,             # kg mol−1 
          "C"   : 0.012011,             # kg mol−1 
          "O"   : 0.015999,             # kg mol−1 
          "N"   : 0.014007,             # kg mol−1 
          "S"   : 0.03206,              # kg mol−1 
          "He"  : 0.0040026,            # kg mol−1 
          "NH3" : 0.017031,             # kg mol−1 
        }

def atm_z(atm, idx):

    # Mean temperature below
    T_mean_down = np.mean( atm.tmp[:idx+1] )

    # # Average temperature weighted by amount of substance?
    # T_mean_down = np.average( atm.tmp[:idx+1], weights=)

    # Use gas phase mean molar mass
    atm.mu[idx] = atm.mu_v[idx]

    # print(atm.grav_z[idx], atm.mu_v[idx], atm.p[idx], atm.p[idx+1])

    # Integration
    dz = - phys.R_gas * T_mean_down * np.log(atm.p[idx+1]/atm.p[idx]) / ( atm.mu[idx] * atm.grav_z[idx] )

    # Next height
    atm.z[idx+1] = atm.z[idx] + dz

    # Next gravity
    atm.grav_z[idx+1] = atm.grav_s * ((atm.planet_radius)**2) / ((atm.planet_radius+atm.z[idx+1])**2)

    # print(idx, T_mean_down, dz, atm.z[idx], atm.z[idx+1], atm.grav_z[idx], atm.grav_z[idx+1])

    return atm

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
        return Tref/(1.-(Tref*phys.R_gas/L_H2O)*math.log(p/pref))
        #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
    if switch == 'CH4':
        Tref = 148.15 # K, arbitrary point (148.15K,esat(148.15K)=9.66bar) on the L/G coexistence curve of methane 
        pref = p_sat('CH4',Tref)
        L_CH4=phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_CH4)*math.log(p/pref))
    if switch == 'CO2':
        Tref = 253. # K, arbitrary point (253K,esat(253K)=20.9bar) on the coexistence curve of CO2 
        pref = p_sat('CO2',Tref)
        L_CO2=phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_CO2)*math.log(p/pref))
    if switch == 'CO':
        Tref = 100. # K, arbitrary point (100K,esat(100K)=4.6bar) on the coexistence curve of CO 
        pref = p_sat('CO',Tref)
        L_CO=phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_CO)*math.log(p/pref))
    if switch == 'N2':
        Tref = 98.15 # K, arbitrary point (98.15K,esat(98.15K)=7.9bar) on the coexistence curve of N2 
        pref = p_sat('N2',Tref)
        L_N2=phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_N2)*math.log(p/pref))
    if switch == 'O2':
        Tref = 123.15 # K, arbitrary point (123.15K,esat(123.15K)=21.9bar) on the coexistence curve of O2 
        pref = p_sat('O2',Tref)
        L_O2=phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_O2)*math.log(p/pref))
    if switch == 'H2':
        Tref = 23.15 # K, arbitrary point (23.15K,esat(23.15K)=1.7bar) on the coexistence curve of H2 
        pref = p_sat('H2',Tref)
        L_H2=phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_H2)*math.log(p/pref))
    if switch == 'He':
        Tref = 4.22 # K, boiling point of He at 1 atm 
        pref = 1e5 # esat('He',Tref) returns 45196 Pa, should return 1e5
        L_He=phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_He)*math.log(p/pref))
    if switch == 'NH3':
        Tref = 273.15 # K, arbitrary point (273.15K,esat(273.15K)=8.6bar) on the coexistence curve of NH3 
        pref = p_sat('NH3',Tref)
        L_NH3=phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_NH3)*math.log(p/pref))
    
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

    # # No latent heat contribution in moist adiabat if below p_sat
    # if P < p_sat(switch, T):
    # # if T > Tdew(switch, psat):
    #     L_heat = 0.

    return L_heat  

# Temperature-dependent molar gas phase heat capacities (J K-1 mol-1)
# https://webbook.nist.gov/chemistry/
def cpv( vol, tmp ):

    # Choose cp functions
    # cp_mode = "constant"    # RTP book
    cp_mode = "T-dependent" # NIST Chemistry WebBook

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
            cp = phys.n2.cp*phys.n2.MolecularWeight*1e-3
        if vol == 'O2':
            cp = phys.o2.cp*phys.o2.MolecularWeight*1e-3
        if vol == 'H2':
            cp = phys.h2.cp*phys.h2.MolecularWeight*1e-3
        if vol == 'He':
            cp = phys.he.cp*phys.he.MolecularWeight*1e-3
        if vol == 'NH3':
            cp = phys.nh3.cp*phys.nh3.MolecularWeight*1e-3   

        return cp # J mol-1 K-1 
'''Adding cp_cond, the heat capacities of the condensates'''
# Temperature-dependent molar condensate heat capacities (J K-1 mol-1)
# https://webbook.nist.gov/chemistry/
def cp_cond( vol, tmp ):


    
    #print(vol)
    # https://webbook.nist.gov/cgi/inchi?ID=C7732185&Mask=2#Thermo-Condensed
    if vol == "H2O":
        # Temperature (K) 298 - 500; 
        A = [-203.6060]
        B = [ 1523.290]
        C = [-3196.413] 
        D = [ 2474.455]
        E = [ 3.855326] 
        F = [-256.5478] 
        G = [-488.7163]
        H = [-285.8304] 
        #if tmp <= 1700:
        cp_idx = 0
        #if tmp > 1700:
        #    cp_idx = 1
        tmp = np.max([tmp, 298]) # Fit validity
        tmp = np.min([tmp,500])
        t = tmp/1000.
        cp = A[cp_idx] + B[cp_idx]*t + C[cp_idx]*t**2. + D[cp_idx]*t**3. + E[cp_idx]/t**2.
    
        return cp # J mol-1 K-1
    # http://www.r744.com/files/pdf_088.pdf
    if vol == "CO2":
        cp_mass = 2048 #J/kg/K <- an uncited value provided in a pdf I found online (see above url); there must be better data on this
        cp_molar = cp_mass*phys.CO2.MolecularWeight/1000 #J/mol/K
        return cp_molar
        
    
    # https://www.osti.gov/etdeweb/servlets/purl/20599211; KAERI Liquid Hydrogen Properties
    if vol == "H2":
        specific_heat_mass_units=14.43877-1.691*tmp + 0.10687*tmp**2-0.00174*tmp**3#J/g/K
        return specific_heat_mass_units*2.02#J/K/mol
    
    # Perkins et al 1991: THE THERMAL CONDUCTIVITY AND HEAT CAPACITY OF FLUID NITROGEN
    if vol == "N2":
        return 60 #J/K/mol; APPROXIMATE VALUE; See Table II or Fig. 3
    
    # https://en.wikipedia.org/wiki/Methane_(data_page)
    if vol == "CH4":
        
        return 52.93 #J/K/mol
    

    # https://www.colby.edu/chemistry/PChem/notes/Ch7Tables.pdf
    if vol == "CO":
        
        return 60 #J/K/mol
    
    # W. F. GIAUQUEA ND H. L. JOHNSTON 1929
    if vol == "O2":
        return 4.184 * 10 #J/K/mol; approximate value, at the higher end of the liquid temperature range

    # No data yet
    if vol == "He":
        return 0. 

    # # SPECIFIC HEAT OF LIQUID AMMONIA by Nathan S. Osborne and Milton S. Van Dusen
    # if vol == "NH3":
    #     mass_specific_heat_capacity = 3.1365 - 0.00057*(tmp+273.15)+16.842/np.sqrt(133-(tmp+273.15)) #joules/gram/K
    #     return mass_specific_heat_capacity*(17.031) #J/k/mol

    # https://www.engineeringtoolbox.com/ammonia-heat-capacity-specific-temperature-pressure-Cp-Cv-d_2016.html#:~:text=At%20ambient%20pressure%20and%20temperature,%5Bcal%2Fg%20K%5D.
    # ! APPLY FIT !
    if vol == "NH3":
        return 4.2239/1000. # J/K/mol
      

satvph2o = phys.satvps_function(phys.H2O)
def slopeRay( logpa, logT ):
    eps     = phys.H2O.MolecularWeight/phys.air.MolecularWeight
    L       = phys.H2O.L_vaporization
    Ra      = phys.air.R
    Rc      = phys.H2O.R
    cpa     = phys.air.cp
    cpc     = phys.H2O.cp
    pa      = math.exp(logpa)
    T       = math.exp(logT)
    qsat    = eps*(satvph2o(T)/pa)
    num     = (1. + (L/(Ra*T))*qsat)*Ra
    den     = cpa + (cpc + (L/(Rc*T) - 1.)*(L/T))*qsat
    return num/den


# Dry adiabat function from p. 91 of RTP book as a comparison
def dry_adiabat( T_surf, P_array, cp_array ):

    T_dry   = np.zeros(len(P_array))
    P_surf  = np.amax(P_array)

    # Check if scalar or cp array
    if not isinstance(cp_array, list):
        cp_array = np.ones(len(P_array))*cp_array

    for idx, prs in enumerate(P_array):

        if cp_array[idx] > 0.:

            T_dry[idx] = T_surf * ( prs / P_surf ) ** ( phys.R_gas / cp_array[idx] )

    return T_dry

# dlnT/dlnP_d, the T-P slope of the non-condensable component of the atmosphere described
# by the general adiabat
def dry_component_slope(lnPd,lnT,atm):
    # T instead lnT
    tmp = math.exp(lnT)

    # Find current atm index
    idx = int(np.amax(atm.ifatm))

    # Sum terms in equation
    num_sum     = 0.
    denom_sum1  = 0. 
    denom_sum2  = 0.
    cp_sum      = 0.
    
    # Calculate sums over volatiles
    for vol in atm.vol_list.keys(): 
        
        # Coefficients
        eta_vol     = atm.x_gas[vol][idx] / atm.xd[idx]
        beta_vol    = L_heat(vol, tmp, atm.p_vol[vol][idx]) / (phys.R_gas * tmp) 

        # Beta terms zero if below saturation vapor pressure
        if atm.p_vol[vol][idx] < p_sat(vol, tmp): beta_vol = 0.

        # Sum in numerator
        num_sum     += eta_vol * beta_vol

        # Sums in denominator
        denom_sum1  += eta_vol * (beta_vol**2.)
        denom_sum2  += -eta_vol * beta_vol
        
        # Adding gaseous cp terms
        cp_sum += eta_vol * cpv(vol, tmp)
        
        # Adding condensate cp terms
        if atm.x_cond[vol][idx] > 0:
            cp_sum += atm.alpha_cloud * atm.x_cond[vol][idx] / atm.xd[idx] * cp_cond(vol, tmp)
            
    # Collect terms
    numerator   = 1. + num_sum
    denominator = (cp_sum / phys.R_gas) + denom_sum1 + denom_sum2

    # dlnT/dlnP_d
    dlnTdlnP_d = numerator / denominator
    
    return dlnTdlnP_d


# dlnT/dlnP slope function 
def moist_slope(lnP, lnT, atm):
    
    # T instead lnT
    tmp = math.exp(lnT)

    # Find current atm index
    idx = int(np.amax(atm.ifatm))

    # Sum terms in equation
    num_sum     = 0.
    denom_sum1  = 0. 
    denom_sum2  = 0. 
    denom_sum3  = 0.

    # Calculate sums over volatiles
    for vol in atm.vol_list.keys(): 
            
        # Coefficients
        eta_vol     = atm.x_gas[vol][idx] / atm.xd[idx]
        beta_vol    = L_heat(vol, tmp, atm.p_vol[vol][idx]) / (phys.R_gas * tmp) 

        # Beta terms zero if below saturation vapor pressure
        if atm.p_vol[vol][idx] < p_sat(vol, tmp): beta_vol = 0.

        # Sum in numerator
        num_sum     += eta_vol * beta_vol

        # Sums in denominator
        denom_sum1  += eta_vol * (beta_vol**2.)
        denom_sum3  += eta_vol
                              
    # Sum 2 in denominator  
    denom_sum2  = num_sum ** 2.

    # Collect terms
    numerator   = 1. + num_sum
    denominator = (atm.cp[idx] / phys.R_gas) + (denom_sum1 + denom_sum2) / (1. + denom_sum3)

    # dlnT/dlnP
    dlnTdlnP = numerator / denominator

    # Moist adiabat slope
    return dlnTdlnP

def condensation( atm, idx, prs_reset):

    # Temperature floor
    tmp = np.amax([atm.tmp[idx], 20.])
    
    # Partial pressures and molar mass: first guess from mixing ratio
    for vol in atm.vol_list.keys():

        # Renormalize mixing ratios to ensure sum == 1
        atm.vol_list[vol]   = atm.vol_list[vol] / sum(atm.vol_list.values())

        # Partial pressures scaled 
        atm.p_vol[vol][idx] = atm.vol_list[vol] * atm.p[idx]

        # Mean gas phase molar mass
        atm.mu_v[idx]       += atm.vol_list[vol] * molar_mass[vol]

    # Account for condensation
    for vol in atm.vol_list.keys():

        # Condensation if p_i > p_sat
        if atm.p_vol[vol][idx] > p_sat(vol, tmp):

            # Set species partial pressure to p_sat
            atm.p_vol[vol][idx]  = p_sat(vol, tmp)

            # Reset mean molar mass, scaled by partial pressure
            mu_v_prev       = atm.mu_v[idx]
            atm.mu_v[idx]   = 0.
            p_vol_sum       = 0.
            for vol1 in atm.vol_list.keys():
                p_vol_sum   += atm.p_vol[vol1][idx]
                atm.mu_v[idx] += molar_mass[vol1] * atm.p_vol[vol1][idx]
            atm.mu_v[idx] /= p_vol_sum

            # Recompute all other partial pressures based on new mean molecular weight
            for vol1 in [ vol1 for vol1 in atm.vol_list.keys() if vol1 != vol and atm.p_vol[vol1][idx] != 0.]:
                
                # print(vol1, atm.p_vol[vol1][idx], end=" ")
                
                # Previous partial pressure of species
                p_vol1_old = atm.p_vol[vol1][idx]

                # New partial pressure due to change in mean molar mass
                p_vol1_new = atm.p_vol[vol1][idx] * ( atm.mu_v[idx] / mu_v_prev )
                
                # Limit to condensation vapor pressure
                atm.p_vol[vol1][idx] = np.min([p_vol1_new, p_sat(vol1, tmp)]) 

                # New total sum of partial pressures
                p_vol_sum += atm.p_vol[vol1][idx] - p_vol1_old

                # P convergence criterion check
                # print(atm.p_vol[vol1][idx], atm.mu_v[idx]/mu_v_prev, p_vol_sum/atm.p[idx], end=" ")


    # Reset mean molar mass
    atm.mu_v[idx]   = 0.
    p_vol_sum       = 0.

    # Update mixing ratios
    for vol in atm.vol_list.keys():

        # Mean molar mass
        p_vol_sum       += atm.p_vol[vol][idx]
        atm.mu_v[idx]   += molar_mass[vol] * atm.p_vol[vol][idx]

        # Condensate phase
        if atm.p_vol[vol][idx] < p_sat(vol, tmp):
            atm.x_cond[vol][idx] = 0.
        else:
            atm.x_cond[vol][idx] = atm.vol_list[vol] - ( atm.p_vol[vol][idx] / atm.p[idx] )

        # atm.x_cond[vol][idx] *= atm.alpha_cloud

        # Add to molar concentration of total condensed phase
        atm.xc[idx]     += atm.x_cond[vol][idx]
        
        # Gas phase molar concentration
        atm.x_gas[vol][idx] = atm.p_vol[vol][idx] / atm.p[idx]
        
        # Add to molar concentration of total gas (dry or moist) phase
        # ! REVISIT ! keeping xd == 0 leads to a bug, why?
        if atm.p_vol[vol][idx] < p_sat(vol, tmp):
            atm.xd[idx]          += atm.x_gas[vol][idx]
        else:
            atm.xv[idx]          += atm.x_gas[vol][idx]
        
        # Mean cp of both gas phase and retained condensates
        atm.cp[idx]    += atm.x_gas[vol][idx] * cpv(vol, atm.tmp[idx]) + atm.x_cond[vol][idx] * cp_cond(vol, atm.tmp[idx]) * atm.alpha_cloud

    # Renormalize cp w/ molar concentration (= cp_hat, Eq. 1, Li+2018)
    atm.cp[idx]  = atm.cp[idx] / ( atm.xd[idx] + atm.xv[idx] )

    # Dry concentration floor
    atm.xd[idx]  = np.amax([atm.xd[idx], 1e-10])
    
    # Correct x values for each volatile; if there's no rain, this changes nothing
    # but if there's rainout, it adjusts them all by the proper fraction
    # NOTE: For some reason this screws up the adiabat a little bit
    atm.xd[idx] *= 1 / ( 1 - ( 1 - atm.alpha_cloud ) * atm.xc[idx] )
    atm.xv[idx] *= 1 / ( 1 - ( 1 - atm.alpha_cloud ) * atm.xc[idx] )
    for vol in atm.vol_list.keys():
        atm.x_cond[vol][idx] *= atm.alpha_cloud / ( 1 - ( 1 - atm.alpha_cloud ) * atm.xc[idx] )
        atm.x_gas[vol][idx]  *= 1 / ( 1 - ( 1 - atm.alpha_cloud) * atm.xc[idx] )
    atm.xc[idx] *= atm.alpha_cloud / ( 1 - ( 1 - atm.alpha_cloud ) * atm.xc[idx] )
    
    return atm

# Builds the generalized moist adiabat from slope(lnP, lnT, atm object)
def general_adiabat( atm ):

    ### Initialization

    # Initialize the tuple solution
    moist_tuple     = [] #[tuple([np.log(atm.ps), atm.ts])] 
   
    # Negative increment to go from ps to ptop < ps       
    step            = -.01

    # Integration counter
    idx             = 0  

    # Surface condensation
    atm             = condensation(atm, idx, prs_reset=False)

    # Create the integrator instance                                              
    int_slope       = integrator(moist_slope, np.log(atm.ps), np.log(atm.ts), step)

    # Update parameters used in the slope function dlntT/dlnP
    int_slope.setParams(atm)
    
    ### Integration of full general adiabat
    while atm.p[idx] > atm.ptop:
        
        # Execute the Runge-Kutta integrator, fill array of tuples
        moist_tuple.append(int_slope.next())

        # Fill next T,P values
        atm.p[idx+1]    = np.exp(int_slope.x)
        atm.tmp[idx+1]  = np.exp(int_slope.y)

        # Calculate next local gravity and height
        atm = atm_z(atm, idx)

        # print("RK4 step, idx:", idx, round(atm.p[idx+1],5), round(atm.tmp[idx+1],5))

        # Set next level to calculate
        idx             += 1
        atm.ifatm[idx]  = idx

        # Calculate condensation at next level
        atm             = condensation(atm, idx, prs_reset=False)

    # Interpolate
    atm = interpolate_atm(atm)

    return atm

# Interpolate and flip pressure, temperature and volatile grids to fixed size
def interpolate_atm(atm):

    # Trim array zeros
    atm_len     = int(np.max(atm.ifatm)+1)
    rest_len    = int(len(atm.p)-atm_len)
    atm.p       = np.flip(np.split(atm.p, [atm_len, rest_len])[0], axis=0)
    atm.tmp     = np.flip(np.split(atm.tmp, [atm_len, rest_len])[0], axis=0)

    # Interpolate staggered nodes
    atm.pl      = np.logspace(np.log10(np.min(atm.p)), np.log10(np.max(atm.p)), atm.nlev_save+1)
    atm.tmpl    = np.interp(atm.pl, atm.p, atm.tmp)

    # Interpolate atmosphere nodes
    prs_itp     = (atm.pl[1:] + atm.pl[:-1]) / 2.
    tmp_itp     = (atm.tmpl[1:] + atm.tmpl[:-1]) / 2.

    # Trim level-dependent quantities
    atm.xd      = np.flip(np.split(atm.xd, [atm_len, rest_len])[0], axis=0)
    atm.xv      = np.flip(np.split(atm.xv, [atm_len, rest_len])[0], axis=0)
    atm.xc      = np.flip(np.split(atm.xc, [atm_len, rest_len])[0], axis=0)
    atm.mrd     = np.flip(np.split(atm.mrd, [atm_len, rest_len])[0], axis=0)
    atm.mrv     = np.flip(np.split(atm.mrv, [atm_len, rest_len])[0], axis=0)
    atm.mrc     = np.flip(np.split(atm.mrc, [atm_len, rest_len])[0], axis=0)
    atm.cp      = np.flip(np.split(atm.cp, [atm_len, rest_len])[0], axis=0)

    # Interpolate level-dependent quantities
    atm.xd      = np.interp(prs_itp, atm.p, atm.xd)
    atm.xv      = np.interp(prs_itp, atm.p, atm.xv)
    atm.xc      = np.interp(prs_itp, atm.p, atm.xc)
    atm.mrd     = np.interp(prs_itp, atm.p, atm.mrd)
    atm.mrv     = np.interp(prs_itp, atm.p, atm.mrv)
    atm.mrc     = np.interp(prs_itp, atm.p, atm.mrc)
    atm.cp      = np.interp(prs_itp, atm.p, atm.cp)

    # Trim & interpolate species-dependent quantities
    for vol in atm.vol_list.keys():
        # Only if volatile is present
        if atm.vol_list[vol] > 1e-10:
            atm.p_vol[vol]   = np.flip(np.split(atm.p_vol[vol], [atm_len, rest_len])[0], axis=0)
            atm.x_gas[vol]   = np.flip(np.split(atm.x_gas[vol], [atm_len, rest_len])[0], axis=0)
            atm.x_cond[vol]  = np.flip(np.split(atm.x_cond[vol], [atm_len, rest_len])[0], axis=0)
            atm.mr_gas[vol]  = np.flip(np.split(atm.mr_gas[vol], [atm_len, rest_len])[0], axis=0)
            atm.mr_cond[vol] = np.flip(np.split(atm.mr_cond[vol], [atm_len, rest_len])[0], axis=0)
    
            atm.p_vol[vol]   = np.interp(prs_itp, atm.p, atm.p_vol[vol])
            atm.x_gas[vol]   = np.interp(prs_itp, atm.p, atm.x_gas[vol]) 
            atm.x_cond[vol]  = np.interp(prs_itp, atm.p, atm.x_cond[vol])
            atm.mr_gas[vol]  = np.interp(prs_itp, atm.p, atm.mr_gas[vol]) 
            atm.mr_cond[vol] = np.interp(prs_itp, atm.p, atm.mr_cond[vol]) 

    # Rewrite atmosphere nodes
    atm.p       = prs_itp
    atm.tmp     = tmp_itp

    return atm

# Plotting
def plot_adiabats(atm):

    # sns.set_style("ticks")
    # sns.despine()

    ls_moist    = 2.5
    ls_dry      = 2.0
    ls_ind      = 1.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))
    # sns.set_style("ticks")
    # sns.despine()
    
    # Rescaling mixing ratios for plotting purposes
    # This should only be uncommented if you remove the rescaling that
    # occurs at the bottom of the condensation() module 
    # atm.xd[:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    # atm.xv[:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    # for vol in atm.vol_list.keys():
    #     if atm.vol_list[vol] > 1e-10:
    #         atm.x_cond[vol][:] *= atm.alpha_cloud / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    #         atm.x_gas[vol][:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    # atm.xc[:] *= atm.alpha_cloud / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    

    # For reference p_sat lines
    T_sat_array    = np.linspace(20,3000,1000) 
    p_partial_sum  = np.zeros(len(atm.tmp))
    
    vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}

    # Individual species
    for vol in vol_list_sorted.keys():

        # Only if volatile is present
        if atm.vol_list[vol] > 1e-10:
    
            # Saturation vapor pressure for given temperature
            Psat_array = [ p_sat(vol, T) for T in T_sat_array ]
            ax1.semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$'+vol_latex[vol], lw=ls_ind, ls=":", color=vol_colors[vol][4])

            # Plot partial pressures
            ax1.semilogy(atm.tmp, atm.p_vol[vol], color=vol_colors[vol][4], lw=ls_ind, ls="-", label=r'$p$'+vol_latex[vol],alpha=0.99)

            # Sum up partial pressures
            p_partial_sum += atm.p_vol[vol]

            # Plot individual molar concentrations
            ax2.semilogy(atm.x_cond[vol],atm.p, color=vol_colors[vol][4], lw=ls_ind, ls="--", label=vol_latex[vol]+" cond.")
            ax2.semilogy(atm.x_gas[vol],atm.p, color=vol_colors[vol][4], lw=ls_ind, ls="-", label=vol_latex[vol]+" gas")
            
    # # Plot sum of partial pressures as check
    ax1.semilogy(atm.tmp, p_partial_sum, color="green", lw=ls_dry, ls="--", label=r'$\sum p^\mathrm{i}$',alpha=0.99)

    # # # Dry adiabat function from RTB book
    # # ! careful: non-T dependent CPs used, can lead to differences with general adiabat !
    # ax1.semilogy( dry_adiabat( atm.ts, atm.pl, atm.cp[-1]), atm.pl , color=vol_colors["black_3"], ls="-.", lw=ls_dry, label=r'Dry adiabat function') # Functional form

    # General moist adiabat
    ax1.semilogy(atm.tmpl, atm.pl, color=vol_colors["black_1"], lw=ls_moist,label="Adiabat",alpha=0.99)

    # Phase molar concentrations
    ax2.semilogy(atm.xd+atm.xv,atm.p, color=vol_colors["black_2"], lw=ls_ind, ls=":", label=r"Gas phase")

    fs_l = 16
    fs_m = 14
    fs_s = 12

    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature, $T$ (K)', fontsize=fs_l)
    ax1.set_ylabel(r'Pressure, $P$ (Pa)', fontsize=fs_l)
    # ax1.set_title('Adiabats & individual Clausius-Clapeyron slopes', fontsize=fs_l)
    ax1.legend(loc=1, ncol=np.min([len(atm.vol_list)+1,2]), fontsize=fs_s)
    ax1.set_xlim([0,np.max(atm.ts)])

    ax2.invert_yaxis()
    # ax2.set_title('Phase & species abundances', fontsize=fs_l)
    ax2.set_xlabel(r'Molar concentration, $X^{\mathrm{i}}_{\mathrm{phase}}$', fontsize=fs_l)
    ax2.set_ylabel(r'Pressure, $P$ (Pa)', fontsize=fs_l)
    ax2.legend(loc=2, ncol=2, fontsize=fs_s)

    ax1.set_ylim(top=atm.ptop)
    ax1.set_ylim(bottom=atm.ps)
    ax2.set_ylim(top=atm.ptop)
    ax2.set_ylim(bottom=atm.ps)

    ax2.set_xscale("log")
    ax2.set_xlim([1e-4, 1.05])
    ax2.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax2.set_xticklabels(["$10^{-4}$", "0.001", "0.01", "0.1", "1"])
    # ax2.set_xlim(right=1.1)

    ax1.tick_params(axis='both', which='major', labelsize=fs_m)
    ax1.tick_params(axis='both', which='minor', labelsize=fs_m)
    ax2.tick_params(axis='both', which='major', labelsize=fs_m)
    ax2.tick_params(axis='both', which='minor', labelsize=fs_m)
    
    ax1.text(0.02, 0.015, 'A', color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l+3, transform=ax1.transAxes)
    ax2.text(0.02, 0.015, 'B', color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l+3, transform=ax2.transAxes)

    #plt.show()

    plt.savefig('./output/general_adiabat.pdf', bbox_inches='tight')
    plt.close(fig)  

    return


####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    # Surface pressure & temperature
    P_surf                  = 260e+5       # Pa
    T_surf                  = 800         # K

    # Volatile molar concentrations: ! must sum to one !
    vol_list = { 
                  "H2O" : .3,    # 300e+5/P_surf --> specific p_surf
                  "CO2" : .3,     # 100e+5/P_surf
                  "H2"  : .3, 
                  "N2"  : .3,     # 1e+5/P_surf
                  "CH4" : .0, 
                  "O2"  : .0, 
                  "CO"  : .0, 
                  "He"  : .0,
                  "NH3" : .0, 
                }

    # Create atmosphere object
    atm                     = atmos(T_surf, P_surf, vol_list)
    
    # Set fraction of condensate retained in column
    atm.alpha_cloud         = 0.0
    
    # Calculate moist adiabat + condensation
    atm                     = general_adiabat(atm)

    # Plot adiabat
    plot_adiabats(atm)
    
