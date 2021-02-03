# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:57:05 2020

@author: sobcr
"""
import numpy as np
import phys
import scipy.interpolate as spint
# Temperature-dependent molar gas phase heat capacities (J K-1 mol-1)
# https://webbook.nist.gov/chemistry/
# Choose cp functions

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
# Thermopedia; 
def cp_cond( vol, tmp, cp_mode='constant'):


    
    #print(vol)
    # 
    
    if vol == "H2O":
        
        temp_array = np.array([273.15,280,285,295,305,315,325,335,345,355,365,373.15])#K
        cp_array = phys.water.MolecularWeight*np.array([4.217,4.198,4.189,4.181,4.178,4.179,4.182,4.186,4.191,4.199,4.209,4.217])#J/K/mol
        cp_interp_func = spint.interp1d(temp_array,cp_array)
        if cp_mode == 'constant':
            cp = phys.water.MolecularWeight * 4.2
        else:
            if tmp < temp_array[0]:
                tmp = temp_array[0]
            elif tmp > temp_array[-1]:
                tmp = temp_array[-1]
            
            cp = cp_interp_func(tmp)
        return cp # J mol-1 K-1
        
        
    # DOI: 10.1615/AtoZ.c.carbon_dioxide
    # https://www.thermopedia.com/content/613/
    if vol == "CO2":
        
        temp_array = np.array([216.6,230,240,250,260,270,280,290,300])#K
        cp_array = phys.CO2.MolecularWeight*np.array([1.84,1.87,1.97,2.1,2.24,2.44,2.81,3.68,8.50])#J/K/mol
        cp_interp_func = spint.interp1d(temp_array,cp_array)
        if cp_mode == 'constant':
            cp = phys.CO2.MolecularWeight * 2.2
        else:
            if tmp < temp_array[0]:
                tmp = temp_array[0]
            elif tmp > temp_array[-1]:
                tmp = temp_array[-1]
            
            cp = cp_interp_func(tmp)
        return cp
        
    
    # https://www.osti.gov/etdeweb/servlets/purl/20599211; KAERI Liquid Hydrogen Properties
    if vol == "H2":
        if cp_mode != 'constant':  
            t = tmp
        else:
            t = 10.
        specific_heat_mass_units=14.43877-1.691*t + 0.10687*t**2-0.00174*t**3#J/g/K
        return specific_heat_mass_units*2.02#J/K/mol
    
    # Perkins et al 1991: THE THERMAL CONDUCTIVITY AND HEAT CAPACITY OF FLUID NITROGEN
    if vol == "N2":
        return 60 #J/K/mol; APPROXIMATE VALUE; See Table II or Fig. 3
    
    # DOI: 10.1615/AtoZ.m.methane
    # https://www.thermopedia.com/content/951/
    if vol == "CH4":
        temp_array = np.array([111.7, 120, 130, 140, 150, 160, 170, 180, 190])#K
        cp_array = phys.CH4.MolecularWeight*np.array([3.48,3.54,3.64,3.80,4.06,4.47,5.23,7.22,92.3])#J/K/mol
        cp_interp_func = spint.interp1d(temp_array,cp_array)
        if cp_mode == 'constant':
            cp = phys.CH4.MolecularWeight * 4.
        else:
            if tmp < temp_array[0]:
                tmp = temp_array[0]
            elif tmp > temp_array[-1]:
                tmp = temp_array[-1]
            
            cp = cp_interp_func(tmp)
        return cp
    

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
    if vol == "NH3":
        return 81.465 # J/K/mol at 298 K
      
