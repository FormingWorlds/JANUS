import numpy as np
import math,phys
import janus.utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from scipy import interpolate
# import seaborn as sns
import copy

# https://webbook.nist.gov/chemistry/
def gas_phase_cp( vol, tmp ):

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
            idx = 0
        if tmp > 1700:
            idx = 1
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
            idx = 0
        if tmp > 1200:
            idx = 1
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
            idx = 0
        if tmp > 1000 and tmp <= 2500:
            idx = 1
        if tmp > 2500:
            idx = 2
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
            idx = 0
        if tmp > 500 and tmp <= 2000:
            idx = 1
        if tmp > 2000:
            idx = 2
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
            idx = 0
        if tmp > 1300:
            idx = 1
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
            idx = 0
        if tmp > 1300:
            idx = 1
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
            idx = 0
        if tmp > 700 and tmp <= 2000:
            idx = 1
        if tmp > 2000:
            idx = 2
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
        idx = 0
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
            idx = 0
        if tmp > 1400:
            idx = 1
        tmp = np.max([tmp, 298]) # Fit validity
        
    t = tmp/1000.
    cp = A[idx] + B[idx]*t + C[idx]*t**2. + D[idx]*t**3. + E[idx]/t**2.

    return cp # J mol-1 K-1


# Vectorize function
gas_phase_cp_vector = np.vectorize(gas_phase_cp)

# Define temperature array
tmp_array = np.linspace(0, 3000, 100)

fig, ax1 = plt.subplots(1, 1, figsize=(7,6))
# sns.set_style("ticks")
# sns.despine()

# Plot temperature vs. cp

for vol in [ "H2O", "CO2", "H2", "N2", "CO", "O2", "He", "NH3" ]:
    ax1.plot(tmp_array,gas_phase_cp_vector(vol, tmp_array), color=ga.vol_colors[vol][4], ls="-", label=ga.vol_latex[vol])

ax1.legend()

ax1.set_xlabel(r'Temperature $T$ (K)')
ax1.set_ylabel(r'c$_\mathrm{p}$ (J mol$^{-1}$ K$^{-1}$)')
# ax4.set_xscale("log")
# ax4.set_yscale("log") 
# ax4.set_xlim(left=0.3, right=100)
# ax4.set_ylim(bottom=1e-20, top=1e5)
# ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])

plt.savefig('./cp_functions.pdf', bbox_inches="tight")
plt.close(fig)

