import copy
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils.SocRadModel import radCompSoc
from modules.dry_adjustment import DryAdj

# Make live plot during time-stepping
def plot_radiative_eqm(new, ref, title):

    unit_bar = 1.0e-5

    fig, ax = plt.subplots(1,1)

    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlabel("Temperature [K]")
    ax.set_title(str(title))

    ax.plot( ref.tmpl, ref.pl*unit_bar, color='grey',  label='$T_{ref}$', ls='--', zorder=2)
    ax.plot( new.tmpl, new.pl*unit_bar, color='black', label='$T_i$', zorder=3)
    
    ax2=ax.twiny()
    ax2.set_xlabel("Heating rate [K/day]", color='tab:red')

    ax2.axvline(x=0, color='tab:red', lw=0.5, alpha=0.8, zorder=0)
    ax2.plot( new.net_heating, new.p*unit_bar, color='tab:red', label='$H_i$', zorder=1)
    ax2_max = np.amax(np.abs(new.net_heating))
    ax2.set_xlim(-1.0*ax2_max, ax2_max)

    fig.legend(loc="lower right")
    fig.savefig("radeqm_iter_PT.png", bbox_inches='tight')  # PNG is faster than PDF
    plt.close()


# Apply heating to atmosphere
def temperature_step(atm, heat, dt, dtmp_clip=1e9, fix_surface=False, smooth=False, dry_adjust=False):

    # Old surface temperature
    tmp_surf_old = atm.tmpl[-1]

    # Calculate temperature change
    dtmp = heat * dt 
    dtmp = np.clip(dtmp, -1.0 * dtmp_clip, dtmp_clip)

    # Apply temperature change
    atm.tmp += dtmp

    # Smooth temperature profile
    if smooth:
        win_len = 5
        pol_odr = 1
        atm.tmp = savgol_filter(atm.tmp, win_len, pol_odr)

    # Dry convective adjustment
    if dry_adjust:
        tmp_before_adj = atm.tmp[:]
        atm = DryAdj(atm)
        adj_changed = np.argwhere(tmp_before_adj != atm.tmp)
        if len(adj_changed) > 0:
            print("    adjusted levels: " + str(adj_changed))

    # Interpolate cell edge temperature values on logp grid
    interp = PchipInterpolator(np.log10(atm.p), atm.tmp, extrapolate=True)
    atm.tmpl = interp(np.log10(atm.pl))

    # Fixed surface temperature?
    if fix_surface:
        atm.tmpl[-1] = tmp_surf_old

    return atm

# Calculate dt array
def calc_stepsize(atm, dt_min=2, dt_max=1e9):

    dtmp_step_frac = 0.5
    dt = dtmp_step_frac * atm.tmp / np.abs(atm.net_heating)
    dt = np.clip(dt, dt_min, dt_max)

    return np.array(dt)


# Find global radiative equilibrium by applying heating rates
def find_radiative_eqm(atm, dirs, rscatter=True, fix_surface=True):

    # Parameters
    steps_max   = 100   # Maximum number of steps
    dtmp_clip   = 10.0  # Maximum change in dT per step [K]
    dtmp_conv   = 5.0   # Change in tmp for convergence [K]
    F_conv      = 0.5   # Flux difference (TOA vs BOA) for convergence [W m-2]
    first_order = False  # Use first-order method
    do_dc_adj   = True  # Do dry convective adjustment

    smooth = True
    F_smooth_stop = 10.0  # Flux error below which to stop smoothing T(p)

    # Variables
    success = False 
    step = 0
    dtmp_dt = np.zeros(np.shape(atm.tmp))
    atm_orig = copy.deepcopy(atm)
    F_loss = np.inf

    # Make isothermal
    atm.tmp[:]  = atm.tmpl[-1]
    atm.tmpl[:] = atm.tmpl[-1]
   
    # Start loop
    while (not success) or (step < steps_max):
        print("    step %d" % step)

        # Store previous atmosphere for reference
        atm_old = copy.deepcopy(atm)
        
        # Get heating rates
        atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
        heat = atm.net_heating

        # Calculate step size to take
        smooth = smooth and (F_loss > F_smooth_stop)
        print("    smooth      = %s" % str(smooth))
        if smooth:
            dt = calc_stepsize(atm, dt_min=2, dt_max=1e9)
        else:
            dt = calc_stepsize(atm, dt_min=0.5, dt_max=10)
        print("    dt_max,med  = %.2f, %.2f days" % (np.amax(dt), np.median(dt)))

        # First order time-stepping method
        if first_order:
            atm_hlf = copy.deepcopy(atm_old)

            # Apply half step
            atm_hlf = temperature_step(atm_hlf, heat, dt * 0.5, fix_surface=fix_surface)

            # Get heating rates at half step
            atm_hlf = radCompSoc(atm_hlf, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
            heat = atm_hlf.net_heating

        # Apply heating rate for full step
        # and smooth temperature profile
        # and do dry convective adjustment
        atm = temperature_step(atm, heat, dt, dtmp_clip=dtmp_clip, fix_surface=fix_surface, smooth=smooth, dry_adjust=do_dc_adj)

        # Calculate statistics
        HR_big_val = np.amax(np.abs(heat))
        dtmp_big_val = np.amax(np.abs(atm.tmp - atm_old.tmp))

        dtmp_dt = np.abs((atm.tmp - atm_old.tmp)/dt)

        F_TOA_net = atm.net_flux[0]
        F_BOA_net = atm.net_flux[-1]
        F_loss = abs(F_TOA_net-F_BOA_net)

        print("    HR_max      = %.3f K/day" % HR_big_val)
        print("    dtmp_max    = %.3f K    " % dtmp_big_val)
        print("    dtmp_dt_max = %.3f K    " % np.amax(dtmp_dt))
        print("    F_TOA_net   = %.1e W m-2" % F_TOA_net)
        print("    F_BOA_net   = %.1e W m-2" % F_BOA_net)

        # Plot
        plt_title = "Step %d \t $|dT_{max}|$ = %.1f K\t $|F_{net}^{loss}|$ = %.1f W m$^{-2}$" % (step, dtmp_big_val, F_loss)
        plot_radiative_eqm(atm, atm_orig, plt_title)

        # Check convergence and prepare new step
        success = np.all(dtmp_dt < dtmp_conv) and (F_loss < F_conv)
        step += 1

        print(" ")


    if not success:
        print("Stopping radiative eqm iterations without success")
    else:
        print("Found global radiative eqm")

    return atm

