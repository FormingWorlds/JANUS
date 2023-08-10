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
def plot_radiative_eqm(atm_hist, ref, title):

    unit_bar = 1.0e-5
    hist_plot = 3

    fig, ax = plt.subplots(1,1)

    # Top axis
    ax2=ax.twiny()
    ax2.set_xlabel("Heating rate [K/day]", color='tab:red')

    ax2.axvline(x=0, color='tab:red', lw=0.5, alpha=0.8, zorder=0)
    ax2.plot( atm_hist[-1].net_heating, atm_hist[-1].p*unit_bar, color='tab:red', label='$H_{%d}(p)$' % int(len(atm_hist)-2), zorder=1)
    ax2_max = np.amax(np.abs(atm_hist[-1].net_heating)) * 1.1
    ax2.set_xlim(-ax2_max, ax2_max)
    ax2.set_xscale("symlog")

    # Bottom axis
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlabel("Temperature [K]")
    ax.set_title(str(title))

    ax.plot( ref.tmpl, ref.pl*unit_bar, color='grey',  label='$T_{ref}(p)$', ls='--', zorder=2)

    plot_count = min(hist_plot, len(atm_hist))
    for i in range(plot_count):
        alpha = (i+1)/plot_count
        idx = len(atm_hist)-plot_count+i
        ax.plot( atm_hist[idx].tmpl, atm_hist[idx].pl*unit_bar, color='black', label='$T_{%d}(p)$' % idx,alpha=alpha)

    # Figure stuff
    fig.legend(loc="lower right")
    fig.savefig("radeqm_iter_PT.png", bbox_inches='tight')  # PNG is faster than PDF
    plt.close()


# Apply heating to atmosphere
def temperature_step(atm, heat, dt, dtmp_clip=1e9, fix_surface=False, smooth=False, adj_steps=0):

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
    if adj_steps > 0:
        tmp_before_adj = copy.deepcopy(atm.tmp)
        for _ in range(adj_steps):
            atm = DryAdj(atm)
        adj_changed = np.count_nonzero(tmp_before_adj - atm.tmp)
        print("    count_adj   = %d" % adj_changed)
        del tmp_before_adj

    # Temperature floor
    atm.tmp = np.clip(atm.tmp,atm.minT, None)

    # Interpolate cell edge temperature values on logp grid
    interp = PchipInterpolator(np.log10(atm.p), atm.tmp, extrapolate=True)
    atm.tmpl = interp(np.log10(atm.pl))

    # Fixed surface temperature?
    if fix_surface:
        atm.tmpl[-1] = tmp_surf_old

    return atm

# Calculate dt array
def calc_stepsize(atm, dt_min=2, dt_max=1e9, dtmp_step_frac=0.8):

    dt = dtmp_step_frac * atm.tmp / np.abs(atm.net_heating)
    dt = np.clip(dt, dt_min, dt_max)

    return np.array(dt)


def find_radiative_eqm(atm, dirs, rscatter=True, surf_state=2, surf_value=0, ini_state=0):
    """Find T(p) satisfying global radiative equilibrium.

    Finds global radiative equilibrium by applying heating rates and dry
    convective adjustment until the difference in net upward flux between
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
            2: Fixed at value in atm object
            3: Fixed at `surf_value`
        ini_state : int
            Initial state of atmosphere.    
            0: isothermal
            1: atm object

    """

    # Parameters
    steps_max   = 250   # Maximum number of steps
    dtmp_clip   = 10.0  # Maximum change in temperature per step [K]
    dtmp_conv   = 10.0   # Change in temperature for convergence [K]
    F_conv      = 0.5   # Flux difference (TOA vs BOA) for convergence [W m-2]
    first_order = False # Use first-order method
    adj_steps   = 3     # Convective adjustment steps (0 for no adjustment)
    fast_start  = True  # Start out fast, with smoothing used to stabilise the model
    F_smooth    = 2.0   # Flux loss below which to stop smoothing
    comp_iter   = 5     # Iterations to look backwards for checking convergence

    # Variables
    success = False 
    step = 1
    atm_hist = []  # Store previous iteration states
    F_loss = np.inf
    smooth = True

    atm_orig = copy.deepcopy(atm)

    # Make isothermal
    if ini_state == 0:
        atm.tmp[:]  = atm.tmpl[-1]
        atm.tmpl[:] = atm.tmpl[-1]

    # Handle surface boundary condition
    match surf_state:
        case 0:
            fix_surface = False 
        case 1:
            fix_surface = False
            print("ERROR: Surface boundary condition %d not implemented")
            exit(1)
        case 2:
            fix_surface = True 
        case 3:
            fix_surface = True
            atm.tmpl[-1] = max(surf_value,atm.minT)

   
    # Start loop
    while (not success) and (step <= steps_max):
        print("    step %d" % step)

        # Store previous atmosphere for reference
        atm_hist.append(copy.deepcopy(atm))
        
        # Get heating rates
        atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
        heat = atm.net_heating

        # Step-size calculation
        smooth = fast_start and smooth and (F_loss > F_smooth)
        print("    smooth      = %s" % str(smooth))
        if smooth:
            dt = calc_stepsize(atm, dt_min=2.0, dt_max=1e9)
        else:
            dt = calc_stepsize(atm, dt_min=0.5, dt_max=1e3, dtmp_step_frac=0.01)
        print("    dt_max,avg  = %.2f, %.2f days" % (np.amax(dt), np.average(dt)))

        # First order time-stepping method (calc HR at half-step)
        if first_order:
            atm_hlf = copy.deepcopy(atm_hist[-1])
            atm_hlf = temperature_step(atm_hlf, heat, dt * 0.5, fix_surface=fix_surface)
            atm_hlf = radCompSoc(atm_hlf, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
            heat = atm_hlf.net_heating

        # Apply heating rate for full step
        # optionally smooth temperature profile
        # optionally do dry convective adjustment
        atm = temperature_step(atm, heat, dt, dtmp_clip=dtmp_clip, fix_surface=fix_surface, smooth=smooth, adj_steps=adj_steps)

        # Calculate statistics
        HR_big = np.amax(np.abs(heat))
        dtmp_dt = np.abs((atm.tmp - atm_hist[-1].tmp)/dt)

        if (step > comp_iter):
            # tmp_compare = np.array([ atm_hist[-comp_iter].tmp, atm_hist[-comp_iter+1].tmp ]).mean(axis=0)
            tmp_compare = atm_hist[-comp_iter].tmp
            dtmp_comp = np.amax(np.abs(atm.tmp - tmp_compare))
        else:
            dtmp_comp = np.inf

        F_TOA_net = atm.net_flux[0]
        F_BOA_net = atm.net_flux[-1]
        F_loss = abs(F_TOA_net-F_BOA_net)

        print("    HR_max      = %.3f K/day" % HR_big)
        print("    dtmp_comp   = %.3f K    " % dtmp_comp)
        print("    dtmp_dt_max = %.3f K    " % np.amax(dtmp_dt))
        print("    F_net^TOA   = %.2e W m-2" % F_TOA_net)
        print("    F_net^BOA   = %.2e W m-2" % F_BOA_net)
        print("    F_net^loss  = %.2f W m-2" % F_loss)

        # Plot
        plt_title = "Step %d     $|dT_{comp}|$ = %.1f K     $|F_{net}^{loss}|$ = %.1f W m$^{-2}$" % (step, dtmp_comp, F_loss)
        plot_radiative_eqm(atm_hist, atm_orig, plt_title)

        # Check convergence
        success = (dtmp_comp < dtmp_conv) and (F_loss < F_conv)
        step += 1

        print(" ")


    if not success:
        print("Stopping radiative equilibrium iterations without success")
    else:
        print("Found global radiative equilibrium")

    return atm

