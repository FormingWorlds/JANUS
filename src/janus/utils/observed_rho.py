import numpy as np

def calc_observed_rho(atm):
    """Calculate the observed bulk density.

    Copied from AGNI.

    Parameters
    ----------
        atm : atmos
            Atmosphere object from atmosphere_column.py

    Returns
    ----------
        rho : float
            Observed bulk density [kg m-3]
    """

    # transspec_r::Float64            # planet radius probed in transmission [m]
    # transspec_m::Float64            # mass [kg] of atmosphere + interior
    # transspec_rho::Float64          # bulk density [kg m-3] implied by r and m

    # arguments
    transspec_p:float = 1e2 # (INPUT) level probed in transmission [Pa]

    # get the observed height
    idx = int(np.argmin(np.abs(atm.p - transspec_p)))
    transspec_r = atm.z[idx] + atm.planet_radius

    # get mass of whole atmosphere, assuming hydrostatic
    transspec_m = np.amax(atm.pl) * 4 * np.pi * atm.planet_radius**2 / atm.grav_s

    # add mass of the interior component
    transspec_m += atm.planet_mass

    # get density of all enclosed by observed layer
    transspec_rho = 3.0 * atm.transspec_m / (4.0 * np.pi * transspec_r**3)

    return transspec_rho

