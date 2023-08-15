# -*- coding: utf-8 -*-
"""
@authors: 
Harrison Nicholls (HN)
"""

import copy, os
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter

from utils.socrates import radCompSoc, CleanOutputDir
from modules.dry_adjustment import DryAdj
from modules.moist_adjustment_H2O import moist_adj

import warnings
warnings.filterwarnings("ignore", message = "All values for SymLogScale")

# Interleave two arrays 
def interleave(a,b):
    # https://stackoverflow.com/questions/5347065/interleaving-two-numpy-arrays-efficiently
    c = np.empty((a.size + b.size,), dtype=a.dtype)
    c[0::2] = a
    c[1::2] = b
    return c

def plot_radiative_eqm(atm_hist, ref, dirs, title, save_as_frame=False):
    """Plot temperature profile and heating rates.

    This plot is used to debug the model as it solves for RC eqm. Can be
    used to make an animation, which is helpful for following model convergence.
    Frames can be combined using ffmpeg:
    `ffmpeg -framerate 5 -i radeqm_monitor_%04d.png -y out.mp4`.

    Parameters
    ----------
        atm_hist : list
            List containing atmosphere objects throughout the model run.
        ref : atmos
            Reference atmosphere object (typically the initial state).
        dirs : dict
            Dictionary of directory paths.
        title : str
            String to be used as the title of the plot
        
        save_as_frame : bool
            Save this image according to frame number, rather than overwriting
            the previous plot. Used for making an animation.
        
    """

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    unit_bar = 1.0e-5
    hist_plot = 3

    fig, (axl, axr) = plt.subplots(1,2, width_ratios=[6, 4])

    # Title
    fig.text(0.5, 0.94, str(title),
        verticalalignment='bottom', horizontalalignment='center',
        transform=fig.transFigure)

    # Top axis
    ax2=axl.twiny()
    ax2.set_xlabel("log(Rad heating [K day-1])", color='tab:red')

    ax2.axvline(x=0, color='tab:red', lw=0.5, alpha=0.8, zorder=0)
    ax2.plot( atm_hist[-1].net_heating, atm_hist[-1].p*unit_bar, color='tab:red', label='$H_{%d}$' % int(len(atm_hist)-1), alpha=0.8)
    
    ax2_max = np.amax(np.abs(atm_hist[-1].net_heating))
    ax2_max = max(ax2_max * 1.5, 1e6)
    ax2.set_xlim(-ax2_max, ax2_max)
    ax2.set_xscale("symlog")
    ax2.set_xticks([])

    # Bottom axis
    axl.set_yscale("log")
    axl.invert_yaxis()
    axl.set_ylabel("Pressure [bar]")
    axl.set_xlabel("Temperature [K]")

    # Right axis
    axr.set_yscale("log")
    axr.invert_yaxis()
    axr.set_xlabel("Net rad flux [W m-2]")
    axr.set_xscale("symlog")
    axr.axvline(x=0, color='black', lw=0.3)

    # Plot data
    axl.plot( ref.tmpl, ref.pl*unit_bar, color='tab:blue',  label='$T_{ref}$', ls='--') # Reference profile
    plot_count = min(hist_plot, len(atm_hist))
    for i in range(plot_count):
        alpha = (i+1)/plot_count
        idx = len(atm_hist)-plot_count+i

        x = interleave(atm_hist[idx].tmpl,atm_hist[idx].tmp) # Include cell edge+centre values
        y = interleave(atm_hist[idx].pl,atm_hist[idx].p) * unit_bar
        axl.plot( x,y, color='black', label='$T_{%d}$' % (idx+1), alpha=alpha)

        axr.plot(atm_hist[idx].net_flux, atm_hist[idx].pl*unit_bar, color='tab:green', alpha=alpha )

    # Legend
    fig.legend(loc="center right")

    # Save figure
    if save_as_frame:
        # Save this file with a unique name so that a video can be made
        # to inspect and debug model convergence. 
        fname = dirs['output']+"/radeqm_monitor_%04d.png" % int(len(atm_hist))
    else:
        # Overwrite last image
        fname = dirs['output']+"/radeqm_monitor.png"

    fig.subplots_adjust(top=0.86)
    fig.savefig(fname, dpi=140)
    plt.close()

def temperature_step(atm, heat, dt, dtmp_clip=1e9, fixed_bottom=True, smooth_width=0, dryadj_steps=0, h2oadj_steps=0):
    """Iterates the atmosphere at each level.

    Includes radiative heating, dry/moist convective adjustment, and smoothing.

    Parameters
    ----------
        atm : atmos
            Atmosphere object
        heat : np.ndarray
            Heating array at each level of the model [K day-1].
        dt : np.ndarray
            Time-step at each level of the model [days].

        dtmp_clip : float
            Maximum absolute change in temperature in this step at each level.
        fixed_bottom : bool
            Fixed bottom edge temperature.
        smooth_width : int
            Convolution window size when performing smoothing. 0 => disabled.
        dryadj_steps : int
            Number of dry convective adjustment steps.
        h2oadj_steps : int
            Number of moist steam convective adjustment steps.
        
    Returns
    ---------- 
        atm : atmos
            Atmosphere object with updated temperature profile
        adj_changed : int
            Number of levels impacted by convective adjustment this iteration.
    """

    # Calculate temperature change
    dtmp = heat * dt 
    dtmp = np.clip(dtmp, -1.0 * dtmp_clip, dtmp_clip)

    # Apply temperature change
    atm.tmp += dtmp

    adj_changed = 0

    # Dry convective adjustment
    if dryadj_steps > 0:
        tmp_before_adj = copy.deepcopy(atm.tmp)
        for _ in range(dryadj_steps):
            atm = DryAdj(atm)
        adj_changed += np.count_nonzero(tmp_before_adj - atm.tmp)
        del tmp_before_adj

    # H2O moist convective adjustment
    if h2oadj_steps > 0:
        tmp_before_adj = copy.deepcopy(atm.tmp)
        atm.tmp += moist_adj(atm, 1.0, nb_convsteps=h2oadj_steps)
        adj_changed += np.count_nonzero(tmp_before_adj - atm.tmp)
        del tmp_before_adj

    # Smooth temperature profile
    if smooth_width > 1:
        atm.tmp = savgol_filter(atm.tmp, smooth_width, 1)

    # Temperature floor
    atm.tmp = np.clip(atm.tmp,atm.minT, None)
    
    # Interpolate cell edge temperature values on logp grid
    bot_old_e = atm.tmpl[-1]
    interp = PchipInterpolator(np.log10(atm.p), atm.tmp, extrapolate=True)
    atm.tmpl = interp(np.log10(atm.pl))

    # Handle boundaries of temperature grid
    atm.tmpl[0] = atm.tmp[0]
    if fixed_bottom:
        atm.tmpl[-1] = bot_old_e
        
    # Second interpolation back to cell-centres, to prevent grid impriting
    # at medium/low resolutions
    interp = PchipInterpolator(np.log10(atm.pl), atm.tmpl, extrapolate=True)
    atm.tmp = interp(np.log10(atm.p))

    return atm, adj_changed

def calc_stepsize(atm, dt_min=1e-5, dt_max=1e7, dtmp_step_frac=1.0):
    """Calculates the time-step at each level for a given accuracy.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

        dt_min : float
            Minimum time-step
        dt_max : float
            Maximum time-step
        dtmp_step_frac : float
            Size of time-step at each level in order to produce the desired
            temperature change, given the current heating rate and temperature.
        
    Returns
    ---------- 
        dt : np.ndarray
            Time-step at each level of the model
    """

    dt_min, dt_max = min(dt_min, dt_max), max(dt_min, dt_max)

    dt = dtmp_step_frac * atm.tmp / np.abs(atm.net_heating)

    dt[-1] *= 0.5 # Reduce dt on bottom layer to avoid numerical problems

    dt = np.array(list(dt))
    dt = np.clip(dt, dt_min, dt_max)

    return dt


def find_rc_eqm(atm, dirs, rscatter=True, 
                surf_state=0, surf_value=350, ini_state=2, 
                gofast=True, dry_adjust=True, h2o_adjust=False, 
                verbose=True, plot=False):
    """Find T(p) satisfying global energy equilibrium.

    Finds global radiative-convective equilibrium by applying heating rates and
    dry convective adjustment until the difference in net upward flux between
    the BOA and TOA is small.

    Parameters
    ----------
        atm : atmos
            Atmosphere object
        dirs : dict
            Directories to save/read things to/from

        rstatter : bool
            Include rayleigh scattering
        surf_state : int
            Temperature boundary condition at the bottom of the model      
            0: Free (although atm.ts is fixed)
            1: Surface balance with turbulent mixing (NOT IMPLEMENTED)
            2: Fixed at `atm.ts`
            3: Fixed at `surf_value`
        surf_value : float 
            Value for BOA temperature when `surf_state==3`
        ini_state : int
            Initial state of atmosphere.    
            2: atm object
            3: isothermal at surface temperature
        gofast : bool
            Speed up convergence by starting out fast, stabilised with smoothing.
        dry_adjust : bool
            Include dry convective adjustment?
        h2o_adjust : bool
            Include pure steam convective adjustment? (Default: False)
        verbose : bool
            Print debug info to stdout?
        plot : bool
            Make plots in output folder?

    Returns
    ----------
        atm : atmos
            Atmosphere object with latest T(p) solution

    """

    # Run parameters
    steps_max    = 100   # Maximum number of steps
    dtmp_gofast  = 35.0  # Change in temperature below which to stop model acceleration
    wait_adj     = 5     # Wait this many steps before introducing convective adjustment
    modprint     = 10    # Frequency to print when verbose==False

    # Convergence criteria
    dtmp_conv    = 10.0   # Maximum rolling change in temperature for convergence (dtmp) [K]
    F_rchng_conv = 0.5   # Maximum relative value of F_loss for convergence [%]
    drel_dt_conv = 5.0    # Maximum rate of relative change in temperature for convergence (dtmp/tmp/dt) [day-1]

    # Variables
    success = False         # Convergence criteria met
    step = 1                # Current step number
    atm_hist = []           # Store previous atmosphere states
    F_loss = 1e99           # Flux loss (TOA vs BOA)
    F_TOA_rad = 1e99        # Net upward TOA radiative flux
    dtmp_comp = np.inf      # Temperature change comparison
    drel_dt = np.inf        # Rate of relative temperature change
    drel_dt_prev  = np.inf  # Previous ^
    flag_prev = False       # Previous iteration is meeting convergence
    step_frac = 1e-3        # Step size fraction relative to absolute temperature
    stopfast = False        # Stopping the 'fast' phase?
    atm_orig = copy.deepcopy(atm)   # Initial atmosphere for plotting

    # Handle initial state
    if ini_state == 3:
        atm.tmp[:]  = atm.ts
        atm.tmpl[:] = atm.ts

    # Handle surface boundary condition
    match surf_state:
        case 0:
            fixed_bottom = False
        case 1:
            fixed_bottom = False
            raise Exception("Atmosphere surface boundary condition %d not implemented")
        case 2:
            fixed_bottom = True
            atm.tmpl[-1] = atm.ts
        case 3:
            fixed_bottom = True
            atm.tmpl[-1] = max(surf_value,atm.minT)
        case _:
            raise ValueError("Invalid surface state for radiative-convective solver")

    # Clean SOCRATES files
    CleanOutputDir(dirs['output'])
    CleanOutputDir( os.getcwd())
   
    # Start loop
    while (not success) and (step <= steps_max):

        # Validate arrays
        if not (np.isfinite(atm.tmp).all() and np.isfinite(atm.tmpl).all()):
            raise Exception("Temperature array contains NaNs")
        if (atm.tmp < 0).any() or (atm.tmpl < 0).any():
            raise Exception("Temperature array contains negative values")

        # Fast initial period
        if gofast and ( (dtmp_comp < dtmp_gofast) or ( step/steps_max > 0.4) ):
            gofast = False 
            stopfast = True
 
        # Step-size calculation and heating rates
        if gofast:
            # Fast phase
            if verbose or (step % modprint == 0): print("    step %d (fast)" % step)
            dtmp_clip = 80.0
            dryadj_steps = 40
            h2oadj_steps = 40
            step_frac = 1e-1
            smooth_window = int( max(0.1*len(atm.tmp),2 )) # 10% of levels

            atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter, rewrite_gas=False, rewrite_cfg=False)
            dt = calc_stepsize(atm, dt_min=1e-2, dtmp_step_frac=step_frac)
            
        else:
            # Slow phase
            if verbose or (step % modprint == 0): print("    step %d" % step)
            dtmp_clip = 10.0
            dryadj_steps = 20
            h2oadj_steps = 20
            step_frac_max = 4e-3
            smooth_window = 0

            if stopfast: # End of 'fast' period (take average of last two iters)
                stopfast = False
                atm.tmp  = np.array([ atm_hist[-1].tmp,  atm_hist[-2].tmp ]).mean(axis=0)
                atm.tmpl = np.array([ atm_hist[-1].tmpl, atm_hist[-2].tmpl]).mean(axis=0)
                dryadj_steps = 0
                h2oadj_steps = 0

            if drel_dt < np.inf:  # Adapt the time-stepping accuracy
                step_frac *= min(max( drel_dt_prev/drel_dt , 0.6 ) , 1.2)
            step_frac = min(step_frac, step_frac_max)
            if verbose: print("    step_frac   = %.2e" % step_frac)

            atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter, rewrite_gas=False, rewrite_cfg=False)
            dt = calc_stepsize(atm, dt_min=1e-5, dt_max=5.0, dtmp_step_frac=step_frac)

        # Cancel convective adjustment if disabled
        if ( dry_adjust and (step < wait_adj)) or not dry_adjust:
            dryadj_steps = 0
        if ( h2o_adjust and (step < wait_adj)) or not h2o_adjust:
            h2oadj_steps = 0

        if verbose: print("    dt_max,med  = %.3f, %.3f days" % (np.amax(dt), np.median(dt)))
        heat = atm.net_heating

        # Apply heating rate for full step
        # optionally smooth temperature profile
        # optionally do dry convective adjustment
        atm, adj_changed = temperature_step(atm, heat, dt, 
                                            dtmp_clip=dtmp_clip, fixed_bottom=fixed_bottom,
                                            smooth_width=smooth_window, 
                                            dryadj_steps=dryadj_steps, h2oadj_steps=h2oadj_steps)

        # Calculate relative rate of change in temperature
        if step > 1:
            drel_dt_prev = drel_dt
            drel_dt = np.amax(np.abs(  ((atm.tmp - atm_hist[-1].tmp)/atm_hist[-1].tmp)/dt  ))

        # Calculate average change in temperature (insensitive to oscillations)
        if step > 3:
            tmp_comp_1 = np.array([ atm.tmp,          atm_hist[-1].tmp ]).mean(axis=0)
            tmp_comp_2 = np.array([ atm_hist[-2].tmp, atm_hist[-3].tmp ]).mean(axis=0)
            dtmp_comp = np.amax(np.abs(tmp_comp_1 - tmp_comp_2))

        # Calculate (the change in) flux balance
        F_TOA_rad_prev = F_TOA_rad
        F_TOA_rad = atm.net_flux[0]
        F_rchng   = abs( (F_TOA_rad - F_TOA_rad_prev) / F_TOA_rad_prev * 100.0)

        F_BOA_rad = atm.net_flux[-1]
        F_OLR_rad = atm.LW_flux_up[0]
        F_loss = abs(F_TOA_rad-F_BOA_rad)

        # Print debug info to stdout
        if verbose:
            # Calculate which level was most impacted
            if step > 1:
                cc_dT = np.amax(np.abs(atm.tmp - atm_hist[-1].tmp))
                ce_dT = np.amax(np.abs(atm.tmpl - atm_hist[-1].tmpl))
                if cc_dT > ce_dT:
                    big_dT_lvl = np.argmax(np.abs(atm.tmp - atm_hist[-1].tmp))
                else:
                    big_dT_lvl = np.argmax(np.abs(atm.tmpl - atm_hist[-1].tmpl)) - 0.5 # i.5 to represent edges
            else:
                big_dT_lvl = -1

            print("    count_adj   = %d layers   " % adj_changed)
            print("    max_change  = %.1f'th lvl " % big_dT_lvl)
            print("    dtmp_comp   = %.3f K      " % dtmp_comp)
            print("    dtmp/tmp/dt = %.3f day-1  " % drel_dt)
            print("    F_rad^OLR   = %.2e W m-2  " % F_OLR_rad)
            print("    F_rad^TOA   = %.2e W m-2  " % F_TOA_rad)
            print("    F_rad^BOA   = %.2e W m-2  " % F_BOA_rad)
            print("    F_rad^loss  = %.2f W m-2  " % F_loss)
            print("    F_chng^TOA  = %.4f %%     " % F_rchng)


        # Store atmosphere for reference
        atm_hist.append(copy.deepcopy(atm))

        # Plot
        if plot:
            plt_title = "Step %d:     $|dT_{comp}|$ = %.1f K     $|F_{rad}^{TOA}|$ = %.1f W m$^{-2}$" % (step, dtmp_comp, F_TOA_rad)
            plot_radiative_eqm(atm_hist, atm_orig, dirs, plt_title, save_as_frame=True)
       
        # Convergence check requires that:
        # - minimal temperature change for two iters
        # - minimal rate of temperature change for two iters
        # - minimal change to net radiative flux at TOA for two iters
        # - solver is not in 'fast' mode, as it is unphysical
        flag_this = (dtmp_comp < dtmp_conv) and (drel_dt < drel_dt_conv) and (F_rchng < F_rchng_conv)
        success   = flag_this and flag_prev
        flag_prev = flag_this and not gofast

        step += 1
        if verbose: print(" ")

    # Print information about the final state
    if not success:
        print("WARNING: Stopping atmosphere iterations without success")
    else:
        print("Convergence criteria met (%d iterations)"%step)

    print("Final radiative fluxes [W m-2]")
    print("    OLR   = %.2e W m-2  " % F_OLR_rad)
    print("    TOA   = %.2e W m-2  " % F_TOA_rad)
    print("    BOA   = %.2e W m-2  " % F_BOA_rad)
    print("    loss  = %.2f W m-2  " % F_loss)

    # Warn user if there's a sign difference in TOA vs BOA fluxes
    if (F_TOA_rad*F_BOA_rad < 0):
        print("WARNING: TOA and BOA radiative fluxes have different signs")

    return atm

