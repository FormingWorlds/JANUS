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

    fig, ax = plt.subplots(1,1)

    # Top axis
    ax2=ax.twiny()
    ax2.set_xlabel("Radiative heating [K day-1]", color='tab:red')

    ax2.axvline(x=0, color='tab:red', lw=0.5, alpha=0.8, zorder=0)
    ax2.plot( atm_hist[-1].net_heating, atm_hist[-1].p*unit_bar, color='tab:red', label='$H_{%d}$' % int(len(atm_hist)-1), alpha=0.8)
    
    ax2_max = np.amax(np.abs(atm_hist[-1].net_heating))
    ax2_max = max(ax2_max * 1.5, 1e6)
    ax2.set_xlim(-ax2_max, ax2_max)
    ax2.set_xscale("symlog")

    # Bottom axis
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlabel("Temperature [K]")
    ax.set_title(str(title))

    ax.plot( ref.tmpl, ref.pl*unit_bar, color='tab:blue',  label='$T_{ref}$', ls='--') # Reference profile
    
    plot_count = min(hist_plot, len(atm_hist))
    for i in range(plot_count):
        alpha = (i+1)/plot_count
        idx = len(atm_hist)-plot_count+i

        x = interleave(atm_hist[idx].tmpl,atm_hist[idx].tmp) # Include cell edge+centre values
        y = interleave(atm_hist[idx].pl,atm_hist[idx].p) * unit_bar
        ax.plot( x,y, color='black', label='$T_{%d}$' % (idx+1), alpha=alpha)

    fig.legend(loc="lower right")

    # Save figure
    if save_as_frame:
        # Save this file with a unique name so that a video can be made
        # to inspect and debug model convergence. 
        fname = dirs['output']+"/radeqm_monitor_%04d.png" % int(len(atm_hist))
    else:
        # Overwrite last image
        fname = dirs['output']+"/radeqm_monitor.png"
    fig.savefig(fname, bbox_inches='tight', dpi=140)  # PNG is faster than PDF
    plt.close()

def temperature_step(atm, heat, dt, dtmp_clip=1e9, fix_surface=False, smooth_width=0, dryadj_steps=0, h2oadj_steps=0):
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
        fix_surface : bool
            Keep surface temperature constant?
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
    interp = PchipInterpolator(np.log10(atm.p), atm.tmp, extrapolate=True)
    atm.tmpl = interp(np.log10(atm.pl))

    # Fixed surface temperature
    if fix_surface:
        atm.tmpl[-1] = atm.ts

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
    dt = np.clip(dt, dt_min, dt_max)

    dt = np.array(list(dt))

    return dt


def find_rc_eqm(atm, dirs, rscatter=True, 
                surf_state=2, surf_value=350, ini_state=2, 
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
            Surface boundary condition.      
            0: Free
            1: Surface balance with turbulent mixing (NOT IMPLEMENTED)
            2: Fixed at `atm.ts`
            3: Fixed at `surf_value`
        surf_value : float 
            Value for surface temperature when `surf_state==3`
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
    second_order = False # Use second order method
    wait_adj     = 5     # Wait this many steps before introducing convective adjustment
    modprint     = 10    # Frequency to print when verbose==False

    # Convergence criteria
    dtmp_conv    = 10.0   # Maximum rolling change in temperature for convergence (dtmp) [K]
    drel_dt_conv = 3.0   # Maximum rate of relative change in temperature for convergence (dtmp/tmp/dt) [day-1]
    F_rloss_conv = 3.0   # Maximum relative value of F_loss for convergence [%]

    # Variables
    success = False         # Convergence criteria met
    step = 1                # Current step number
    atm_hist = []           # Store previous atmosphere states
    F_loss = 1e99           # Flux loss (TOA vs BOA)
    F_loss_prev = 1e99      # Previous ^
    dtmp_comp = np.inf      # Temperature change comparison
    drel_dt = np.inf        # Rate of relative temperature change
    drel_dt_prev  = np.inf  # Previous ^
    step_frac = 0.01        # Step size scaling relative to temperature
    flag_prev = False       # Previous iteration is meeting convergence
    stopfast = False        # Stopping the 'fast' phase?
    atm_orig = copy.deepcopy(atm)   # Initial atmosphere for plotting

    # Handle initial state
    if ini_state == 3:
        atm.tmp[:]  = atm.ts
        atm.tmpl[:] = atm.ts

    # Handle surface boundary condition
    match surf_state:
        case 0:
            fix_surface = False 
        case 1:
            fix_surface = False
            raise Exception("Atmosphere surface boundary condition %d not implemented")
        case 2:
            fix_surface = True 
            atm.tmpl[-1] = atm.ts
        case 3:
            fix_surface = True
            atm.tmpl[-1] = max(surf_value,atm.minT)
            atm.ts = atm.tmpl[-1]
        case _:
            raise ValueError("Invalid surface state for radiative-convective solver")

    # Clean SOCRATES files
    CleanOutputDir(dirs['output'])
    CleanOutputDir( os.getcwd())
   
    # Start loop
    while (not success) and (step <= steps_max):

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
            smooth_window = 10

            atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter, rewrite_gas=False, rewrite_cfg=False)
            dt = calc_stepsize(atm, dt_min=1e-2, dtmp_step_frac=0.1)
            
        else:
            # Slow phase
            if verbose or (step % modprint == 0): print("    step %d" % step)
            dtmp_clip = 10.0
            dryadj_steps = 20
            h2oadj_steps = 20
            smooth_window = 0

            if drel_dt_prev < np.inf:  # Adapt the time-stepping accuracy
                step_frac *= min(max( (drel_dt_prev/drel_dt) , 0.6 ) , 1.2)
            step_frac = min(step_frac, 5e-3)
            if verbose: print("    step_frac   = %.2e" % step_frac)
            
            if stopfast: # End of 'fast' period (take average of last two iters)
                stopfast = False
                atm.tmp  = np.array([ atm_hist[-1].tmp,  atm_hist[-2].tmp ]).mean(axis=0)
                atm.tmpl = np.array([ atm_hist[-1].tmpl, atm_hist[-2].tmpl]).mean(axis=0)
                dryadj_steps = 0
                h2oadj_steps = 0
                smooth_window = 10

            atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter, rewrite_gas=False, rewrite_cfg=False)
            dt = calc_stepsize(atm, dt_max=50, dtmp_step_frac=step_frac)

        # Cancel convective adjustment if disabled
        if ( dry_adjust and (step < wait_adj)) or not dry_adjust:
            dryadj_steps = 0
        if ( h2o_adjust and (step < wait_adj)) or not h2o_adjust:
            h2oadj_steps = 0

        if verbose: print("    dt_max,med  = %.3f, %.3f days" % (np.amax(dt), np.median(dt)))
        heat = atm.net_heating

        # Optional second order time-stepping method (calculate HR at half-step)
        if second_order:
            atm_hlf = copy.deepcopy(atm_hist[-1])
            atm_hlf = temperature_step(atm_hlf, heat, dt * 0.5, dtmp_clip=dtmp_clip, fix_surface=fix_surface)
            if surf_state == 0:
                atm_hlf.ts = atm_hlf.tmpl[-1]
            atm_hlf = radCompSoc(atm_hlf, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
            heat = atm_hlf.net_heating

        # Apply heating rate for full step
        # optionally smooth temperature profile
        # optionally do dry convective adjustment
        atm, adj_changed = temperature_step(atm, heat, dt, dtmp_clip=dtmp_clip, 
                                            fix_surface=fix_surface, smooth_width=smooth_window, 
                                            dryadj_steps=dryadj_steps, h2oadj_steps=h2oadj_steps)
        if surf_state == 0:
            atm.ts = atm.tmpl[-1]

        # Calculate relative rate of change in temperature
        if step > 2:
            drel_dt_prev = drel_dt
            drel_dt = np.amax(np.abs(  ((atm.tmp - atm_hist[-1].tmp)/atm_hist[-1].tmp)/dt  ))

        # Calculate average change in temperature (insensitive to oscillations)
        if step > 3:
            tmp_comp_1 = np.array([ atm.tmp,          atm_hist[-1].tmp ]).mean(axis=0)
            tmp_comp_2 = np.array([ atm_hist[-2].tmp, atm_hist[-3].tmp ]).mean(axis=0)
            dtmp_comp = np.amax(np.abs(tmp_comp_1 - tmp_comp_2))

        # Calculate (the change in) flux balance
        F_TOA_rad = atm.net_flux[0]
        F_BOA_rad = atm.net_flux[-1]
        F_OLR_rad = atm.LW_flux_up[0]
        F_loss_prev = F_loss
        F_loss = abs(F_TOA_rad-F_BOA_rad)
        F_rloss = abs(F_loss - F_loss_prev) / F_loss_prev * 100.0

        if verbose:
            print("    count_adj   = %d layers   " % adj_changed)
            print("    dtmp_comp   = %.3f K      " % dtmp_comp)
            print("    dtmp/tmp/dt = %.3f day-1  " % drel_dt)
            print("    F_rad^TOA   = %.2e W m-2  " % F_TOA_rad)
            print("    F_rad^BOA   = %.2e W m-2  " % F_BOA_rad)
            print("    F_rad^OLR   = %.2e W m-2  " % F_OLR_rad)
            print("    F_rad^loss  = %.2f W m-2  " % F_loss)
            print("    F_rad^rloss = %.3f %%     " % F_rloss)

        # Store atmosphere for reference
        atm_hist.append(copy.deepcopy(atm))

        # Plot
        if plot:
            plt_title = "Step %d:     $|dT_{comp}|$ = %.1f K     $|F_{rad}^{loss}|$ = %.1f W m$^{-2}$" % (step, dtmp_comp, F_loss)
            plot_radiative_eqm(atm_hist, atm_orig, dirs, plt_title, save_as_frame=True)
       
        # Convergence check requires that:
        # - minimal temperature change for two iters
        # - minimal rate of temperature change for two iters
        # - minimal change to radiative flux loss for two iters
        # - solver is not in 'fast' mode, as it is unphysical
        flag_this = (dtmp_comp < dtmp_conv) and (drel_dt < drel_dt_conv) and (F_rloss < F_rloss_conv)
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

    return atm

