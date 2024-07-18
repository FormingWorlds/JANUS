import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from janus.utils.atmosphere_column import atmos
import janus.utils.phys as phys

def plot_emission(atm:atmos, filename:str='output/toa_emission.pdf', 
                  level_idx:int=0, planck_surface:bool=True, 
                  incoming_solar:bool=True, show_bands:bool=False):

    # Configure
    planck_samps = 1000
    lw = 0.5
    xmin = 0.1 # nm
    xmax = 3e5   # nm

    # Get data
    x = np.zeros(atm.nbands)
    y = np.zeros(atm.nbands)

    for i in range(atm.nbands):

        # Wavelength
        x[i] += atm.band_centres[i]

        # Flux 
        y[i] += atm.LW_spectral_flux_up[i,level_idx]  # W m-2
        y[i] += atm.SW_spectral_flux_up[i,level_idx]  # W m-2
        y[i] /= atm.band_widths[i]  # W m-2 nm-1
        y[i] *= 1e3  # erg s-1 cm-2 nm-1


    # Make plot 
    fig,ax = plt.subplots(1,1, figsize=(8,4))

    # Plot incoming stellar spectrum
    if incoming_solar:
        ysolar = atm.SW_spectral_flux_down.T[0] # W m-2 
        ysolar /= atm.band_widths 
        ysolar *= 1e3 
        ax.plot(x, ysolar, color="green", lw=lw, label="TOA stellar flux")

    # Plot band edges
    if show_bands:
        for i in range(atm.nbands):
            c = 'pink'
            l = ""
            if i == 0:
                l = "Band edges"
            ax.axvline(x=atm.band_edges[i]  , color=c, lw=0.2, label=l)
            ax.axvline(x=atm.band_edges[i+1], color=c, lw=0.2)
    
    # Plot surface emission
    if planck_surface:
        planck_tmp = atm.ts
        xp = np.logspace(np.log10(np.amin(x))  ,  np.log10(np.amax(x))  ,  planck_samps)
        yp = np.zeros(np.shape(xp))
        for i in range(planck_samps):
            wl = xp[i] * 1.0e-9 # metres

            # Calculate planck function value [W m-2 sr-1 m-1]
            # http://spiff.rit.edu/classes/phys317/lectures/planck.html
            yp[i] = 2.0 * phys.h * phys.c**2.0 / wl**5.0   *   1.0 / ( np.exp(phys.h * phys.c / (wl * phys.k * planck_tmp)) - 1.0)

            # Integrate solid angle (hemisphere), scale by albedo, convert units
            yp[i] = yp[i] * np.pi * 1.0e-9 # [W m-2 nm-1]
            yp[i] = yp[i] * (1.0 - atm.albedo_s)
            yp[i] = yp[i] * 1000.0 # [erg s-1 cm-2 nm-1]
        ax.plot(xp,yp, color="dodgerblue", lw=lw, label="Surface")
    
    # Plot spectrum
    ax.plot(x,y, color='k', lw=lw, label="Emission ($i=%d$)"%level_idx)

    # Adjust
    ax.set_xlabel("Wavelength [nm]")
    ax.set_xlim(max(xmin, np.amin(x)), min(xmax, np.amax(x)))
    ax.set_xscale("log")

    ax.set_ylabel("Flux density [erg s$^{-1}$ cm$^{-2}$ nm$^{-1}$]")
    ax.set_yscale("log")

    ax.legend(loc="lower center")

    wlwn = lambda x: 1.0e7/(x-1.0e-20)  # wavelength <-> wavenumber, with small offset to prevent divergence
    secax = ax.secondary_xaxis('top', functions=(wlwn, wlwn))
    secax.set_xlabel('Wavenumber [cm$^{-1}$]')

    # Save fig 
    fig.savefig(filename, bbox_inches='tight', dpi=190)
    plt.close()

