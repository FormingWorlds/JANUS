import numpy as np

# Triple point of water
T_tp = 273.16
# Critical temperature of water
T_c = 647.096
# Number of points in table
N = 50

# Temperature array on which data is defined
T = np.array([273.16   , 280.79134694, 288.42269388, 296.05404082,
       303.68538776, 311.31673469, 318.94808163, 326.57942857,
       334.21077551, 341.84212245, 349.47346939, 357.10481633,
       364.73616327, 372.3675102 , 379.99885714, 387.63020408,
       395.26155102, 402.89289796, 410.5242449 , 418.15559184,
       425.78693878, 433.41828571, 441.04963265, 448.68097959,
       456.31232653, 463.94367347, 471.57502041, 479.20636735,
       486.83771429, 494.46906122, 502.10040816, 509.7317551 ,
       517.36310204, 524.99444898, 532.62579592, 540.25714286,
       547.8884898 , 555.51983673, 563.15118367, 570.78253061,
       578.41387755, 586.04522449, 593.67657143, 601.30791837,
       608.93926531, 616.57061224, 624.20195918, 631.83330612,
              639.46465306, T_c])

# Latent heat of vaporization [J/kg]
L_vap = np.array([2500537.992720097, 2482839.9845154616,2464924.549801353,
       2446870.4063889873,2428716.9993550396,2410476.4281226317,
       2392141.8576981686,2373693.4045157116,2355102.225994891,
       2336333.3440693025,2317347.58344357,2298102.8941141665,
       2278555.2457944946,2258659.222272561,2238368.4009194765,
       2217635.5723345308,2196412.8342096675,2174651.5794327096,
       2152302.389309872,2129314.8371038437,2105637.2037150585,
       2081216.1053961227,2055996.0321805072,2029918.7946571708,
       2002922.8753449027,1974942.6787642944,1945907.670911796,
       1915741.3936820417,1884360.3321936703,1851672.602035535,
       1817576.4078452894,1781958.2023893306,1744690.443404617,
       1705628.7991290821,1664608.5850239818,1621440.1110193927,
       1575902.4591775925,1527734.9582067847,1476625.2050181734,
       1422191.7739563962, 1363958.4911911818,1301314.785717382,
       1233451.918479445,1159254.7924653592,1077105.2789987705,
       984489.844927971,877106.1104059853,746366.06771214,
                  569229.6959539498, 0.0])

# Saturation vapour pressure [Pa]
psat=np.array([6.11657070e+02, 1.04700395e+03, 1.73586788e+03, 2.79466788e+03,
       4.37920771e+03, 6.69301793e+03, 9.99615270e+03, 1.46142123e+04,
       2.09473490e+04, 2.94790161e+04, 4.07842377e+04, 5.55372020e+04,
       7.45180209e+04, 9.86185378e+04, 1.28847110e+05, 1.66332336e+05,
       2.12325737e+05, 2.68203431e+05, 3.35466882e+05, 4.15742812e+05,
       5.10782391e+05, 6.22459837e+05, 7.52770551e+05, 9.03828931e+05,
       1.07786602e+06, 1.27722709e+06, 1.50436941e+06, 1.76186025e+06,
       2.05237530e+06, 2.37869784e+06, 2.74371868e+06, 3.15043722e+06,
       3.60196390e+06, 4.10152446e+06, 4.65246639e+06, 5.25826816e+06,
       5.92255199e+06, 6.64910116e+06, 7.44188315e+06, 8.30508042e+06,
       9.24313145e+06, 1.02607859e+07, 1.13631798e+07, 1.25559402e+07,
       1.38453372e+07, 1.52385141e+07, 1.67438702e+07, 1.83717818e+07, 
       2.01363422e+07, 2.20640000e+07])

  #dln(psat)/dln(T) -- equivalent to L/R_v/T when constant L
phase_grad = np.array([19.84498763, 19.17337147, 
       18.53757655, 17.93520481, 17.36408297,
       16.82223398, 16.30785266, 15.819285  , 15.35501047, 14.91362679,
       14.49383689, 14.09443772, 13.71431049, 13.35241236, 13.00776922,
       12.6794694 , 12.36665832, 12.06853382, 11.78434214, 11.51337446,
       11.254964  , 11.00848346, 10.77334297, 10.54898836, 10.33489974,
       10.13059049,  9.93560646,  9.74952556,  9.57195763,  9.40254465,
        9.24096135,  9.08691622,  8.940153  ,  8.80045284,  8.66763718,
        8.54157171,  8.42217156,  8.30940852,  8.20332077,  8.1040267 ,
        8.01174465,  7.92682246,  7.84978332,  7.78140111,  7.72283313,
                       7.67587698,  7.64353758,  7.63155818,  7.65437086, 7.85951783])

tables = {'T':T,
          'L_vap':L_vap,
          'psat':psat,
          'phase_grad':phase_grad
          }

def lookup(table, T_val):
    """
    HII 6/23
    --------------------------------------------------------------------
    Lookup a single value for one of the thermodynamic properties of water
    Uses fact that T values are linearly spaced. Should work for either
    single T_val or numpy array of T_val
    ---------------------------------------------------------------------
    Input:
    table, type: Str -- name of one of the lookup tables
    T_val, type: float -- temperature(s) at which to interpolate result
    ----------------------------------------------------------------------
    Returns:
    Value of thermodynamic property of water interpolated to temp T from 
    table 'table'
    """
    data = tables[table]
    i =np.floor((T_val - T[0])/(T[-1] - T[0]) * (N-1) ).astype(int)

    print('hello', np.amax(i), np.any(i>N-1))
    
    if np.any(i<0) or np.any(i>N-2):
        raise Exception('Error: T_val given that is out of range')
        
    return data[i] + (T_val - T[i])/(T[i+1] - T[i]) * (data[i+1] - data[i])

    
