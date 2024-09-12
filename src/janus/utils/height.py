import numpy as np
import janus.utils.phys as phys

import logging
log = logging.getLogger("fwl."+__name__)


def gravity( m, r ):
    g = phys.G*m/r**2
    return g

def integrate_heights(atm, m_planet, r_planet):

    z       = np.zeros(len(atm.p))
    grav_s          = gravity( m_planet, r_planet )

    # Reverse arrays to go from high to low pressure
    atm.p   = atm.p[::-1]
    atm.tmp = atm.tmp[::-1]
    for vol in atm.vol_list.keys():
        atm.x_gas[vol] = atm.x_gas[vol][::-1]

    height_error = False
    for n in range(0, len(z)-1):

        # Gravity with height
        grav_z = grav_s * ((r_planet)**2) / ((r_planet + z[n])**2)

        # Mean molar mass depending on mixing ratio
        mean_molar_mass = 0
        for vol in atm.vol_list.keys():
            mean_molar_mass += phys.molar_mass[vol]*atm.x_gas[vol][n]

        # Use hydrostatic equation to get height difference
        dz = phys.R_gas * atm.tmp[n] / (mean_molar_mass * grav_z * atm.p[n]) * (atm.p[n] - atm.p[n+1])

        # Next height
        z[n+1] = z[n] + dz

        # Check if heights are very large.
        # This implies that the hydrostatic/gravity integration failed.
        if (z[n+1] > 1.0e8) or (dz > 1e8):
            height_error = True
            log.error("Hydrostatic integration blew up. Setting dummy values for height")
            break

    # Set dummy values
    if height_error:
        z = np.linspace(0.0, 1000.0, len(atm.p))

    # Reverse arrays again back to normal
    atm.p   = atm.p[::-1]
    atm.tmp = atm.tmp[::-1]
    for vol in atm.vol_list.keys():
        atm.x_gas[vol] = atm.x_gas[vol][::-1]
    z = z[::-1]

    # Set cell edge values
    zl = np.zeros(len(z)+1)
    for i in range(1,len(z)):
        zl[i] = 0.5 * (z[i-1] + z[i])
    zl[0] = 2*z[0] - zl[1] # estimate TOA height

    return z, zl, height_error
