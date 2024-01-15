import numpy as np
import utils.phys as phys

def gravity( m, r ):
    g = phys.G*m/r**2
    return g

def AtmosphericHeight(atm, m_planet, r_planet):

    z_profile       = np.zeros(len(atm.p))
    P_s             = np.max(atm.p)
    grav_s          = gravity( m_planet, r_planet )

    # Reverse arrays to go from high to low pressure
    atm.p   = atm.p[::-1]
    atm.tmp = atm.tmp[::-1]
    for vol in atm.vol_list.keys():
        atm.x_gas[vol] = atm.x_gas[vol][::-1]

    # print(atm.p)

    for n in range(0, len(z_profile)-1):

        # Gravity with height
        grav_z = grav_s * ((r_planet)**2) / ((r_planet + z_profile[n])**2)

        # print(r_planet, grav_s, grav_z, z_profile[n])

        # Mean molar mass depending on mixing ratio
        mean_molar_mass = 0
        for vol in atm.vol_list.keys():
            mean_molar_mass += phys.molar_mass[vol]*atm.x_gas[vol][n]

        # Temperature below present height
        # T_mean_below    = np.mean(atm.tmp[n:])

        # # Direction calculation
        # z_profile[n] = - R_gas * T_mean_below * np.log(atm.p[n]/P_s) / ( mean_molar_mass * grav_s )

        # Integration
        # dz = - phys.R_gas * T_mean_below * np.log(atm.p[n+1]/atm.p[n]) / (mean_molar_mass*grav_z)
        dz = phys.R_gas * atm.tmp[n] / (mean_molar_mass * grav_z * atm.p[n]) * (atm.p[n] - atm.p[n+1]) 
        
        # Next height
        z_profile[n+1] = z_profile[n] + dz

    # Reverse arrays again back to normal
    atm.p   = atm.p[::-1]
    atm.tmp = atm.tmp[::-1]
    for vol in atm.vol_list.keys():
        atm.x_gas[vol] = atm.x_gas[vol][::-1]
    z_profile = z_profile[::-1]

    return z_profile
