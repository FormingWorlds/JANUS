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
def plot_radiative_eqm(atm_hist, ref, dirs, title, save_frames=False):

    unit_bar = 1.0e-5
    hist_plot = 3

    fig, ax = plt.subplots(1,1)

    # Top axis
    ax2=ax.twiny()
    ax2.set_xlabel("Radiative heating rate [K/day]", color='tab:red')

    ax2.axvline(x=0, color='tab:red', lw=0.5, alpha=0.8, zorder=0)
    ax2.plot( atm_hist[-1].net_heating, atm_hist[-1].p*unit_bar, color='tab:red', label='$H_{%d}$' % int(len(atm_hist)), alpha=0.8)
    
    ax2_max = int(np.amax(np.abs(atm_hist[-1].net_heating)) * 1.5) + 1.0
    ax2.set_xlim(-ax2_max, ax2_max)
    ax2.set_xscale("symlog")

    # Bottom axis
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlabel("Temperature [K]")
    ax.set_title(str(title))

    ax.plot( ref.tmpl, ref.pl*unit_bar, color='tab:blue',  label='$T_{ref}$', ls='--')

    plot_count = min(hist_plot, len(atm_hist))
    for i in range(plot_count):
        alpha = (i+1)/plot_count
        idx = len(atm_hist)-plot_count+i
        ax.plot( atm_hist[idx].tmpl, atm_hist[idx].pl*unit_bar, color='black', label='$T_{%d}$' % (idx+1), alpha=alpha)

    # Figure stuff
    fig.legend(loc="lower right")

    if save_frames:
        fname = dirs['output']+"/radeqm_monitor.%d.png" % int(len(atm_hist))
    else:
        fname = dirs['output']+"/radeqm_monitor.png"
    fig.savefig(fname, bbox_inches='tight', dpi=140)  # PNG is faster than PDF
    plt.close()

# Apply heating to atmosphere
def temperature_step(atm, heat, dt, dtmp_clip=1e9, fix_surface=False, smooth_width=0, adj_steps=0):

    # Calculate temperature change
    dtmp = heat * dt 
    dtmp = np.clip(dtmp, -1.0 * dtmp_clip, dtmp_clip)

    # Apply temperature change
    atm.tmp += dtmp

    # Dry convective adjustment
    if adj_steps > 0:
        tmp_before_adj = copy.deepcopy(atm.tmp)
        for _ in range(adj_steps):
            atm = DryAdj(atm)
        adj_changed = np.count_nonzero(tmp_before_adj - atm.tmp)
        print("    count_adj   = %d layers" % adj_changed)
        del tmp_before_adj

    # Smooth temperature profile
    if smooth_width > 1:
        atm.tmp = savgol_filter(atm.tmp, smooth_width, 1)

    # Temperature floor
    atm.tmp = np.clip(atm.tmp,atm.minT, None)

    # Interpolate cell edge temperature values on logp grid
    interp = PchipInterpolator(np.log10(atm.p), atm.tmp, extrapolate=True)
    atm.tmpl = interp(np.log10(atm.pl))

    # Fixed surface temperature?
    if fix_surface:
        atm.tmpl[-1] = atm.ts

    return atm

# Calculate dt array
def calc_stepsize(atm, dt_min=1e-5, dt_max=1e9, dtmp_step_frac=1.0):

    dt = dtmp_step_frac * atm.tmp / np.abs(atm.net_heating)
    dt = np.clip(dt, dt_min, dt_max)

    return np.array(dt)


def find_radiative_eqm(atm, dirs, rscatter=True, surf_state=2, surf_value=350, ini_state=2, gofast=True):
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
            2: Fixed at `atm.ts`
            3: Fixed at `surf_value`
        surf_value : float 
            Value for surface temperature when `surf_state==3`
        ini_state : int
            Initial state of atmosphere.    
            2: atm object
            3: isothermal at surface temperature
        gofast : bool
            Speed up convergence by starting out fast with smoothing to stabilise

    """

    # Parameters
    steps_max    = 100   # Maximum number of steps
    dtmp_conv    = 10.0  # Maximum change in temperature for convergence [K]
    drel_dt_conv = 5.0   # Maximum rate of relative change in temperature for convergence (dtmp/tmp/dt) [days-1]
    adj_steps    = 40    # Convective adjustment steps (0 for no adjustment)
    dtmp_gofast  = 20.0  # Change in temperature below which to stop model acceleration

    # Variables
    success = False 
    step = 1
    atm_hist = []  # Store previous iteration states
    F_loss = np.inf
    dtmp_comp = np.inf
    drel_dt = np.inf
    drel_dt_prev  = np.inf
    step_frac = 0.01
    dtmp_dt = np.inf
    flag_previous = False
    atm_orig = copy.deepcopy(atm)

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
            print("ERROR: Surface boundary condition %d not implemented")
            exit(1)
        case 2:
            fix_surface = True 
            atm.tmpl[-1] = atm.ts
        case 3:
            fix_surface = True
            atm.tmpl[-1] = max(surf_value,atm.minT)
   
    # Start loop
    while (not success) and (step <= steps_max):
        print("    step %d" % step)

        # Get heating rates
        atm = radCompSoc(atm, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
        heat = atm.net_heating

        # Step-size calculation
        gofast = gofast and (dtmp_comp > dtmp_gofast)
        if gofast:
            print("    gofast      = True")
            first_order = False
            dtmp_clip = 80.0
            smooth_window = 12
            dt = calc_stepsize(atm, dt_min=2.0, dtmp_step_frac=0.4)
            
        else:
            first_order = True
            dtmp_clip   = 30.0
            smooth_window = 0

            if drel_dt_prev < np.inf:
                step_frac *= min(max( (drel_dt_prev/drel_dt) , 0.6 ) , 1.2)
            step_frac = min(step_frac, 4e-3)
            print("    step_frac   = %.2e" % step_frac)
            dt = calc_stepsize(atm, dt_max=0.4, dtmp_step_frac=step_frac)

        
        print("    dt_max,med  = %.3f, %.3f days" % (np.amax(dt), np.median(dt)))

        # First order time-stepping method (calc HR at half-step)
        if first_order:
            atm_hlf = copy.deepcopy(atm_hist[-1])
            atm_hlf = temperature_step(atm_hlf, heat, dt * 0.5, dtmp_clip=dtmp_clip, fix_surface=fix_surface, adj_steps=adj_steps)
            atm_hlf = radCompSoc(atm_hlf, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
            heat = atm_hlf.net_heating

        # Apply heating rate for full step
        # optionally smooth temperature profile
        # optionally do dry convective adjustment
        atm = temperature_step(atm, heat, dt, dtmp_clip=dtmp_clip, fix_surface=fix_surface, smooth_width=smooth_window, adj_steps=adj_steps)

        # Calculate statistics
        HR_max = np.amax(np.abs(heat))
        if step > 2:
            dtmp_dt = np.amax(np.abs((atm.tmp - atm_hist[-1].tmp)/dt))
            drel_dt_prev = drel_dt
            drel_dt = np.amax(np.abs((atm.tmp - atm_hist[-1].tmp)/atm_hist[-1].tmp/dt))

        if (step > 5):
            tmp_comp_1 = np.array([ atm.tmp,          atm_hist[-1].tmp ]).mean(axis=0)
            tmp_comp_2 = np.array([ atm_hist[-2].tmp, atm_hist[-3].tmp ]).mean(axis=0)
            dtmp_comp = np.amax(np.abs(tmp_comp_1 - tmp_comp_2))

        F_TOA_rad = atm.net_flux[0]
        F_BOA_rad = atm.net_flux[-1]
        F_loss = abs(F_TOA_rad-F_BOA_rad)

        print("    dtmp_comp   = %.3f K      " % dtmp_comp)
        print("    HR_max      = %.3f K day-1" % HR_max)
        print("    dtmp/dt     = %.3f K day-1" % dtmp_dt)
        print("    dtmp/tmp/dt = %.3f day-1  " % drel_dt)
        print("    F_rad^TOA   = %.2e W m-2  " % F_TOA_rad)
        print("    F_rad^BOA   = %.2e W m-2  " % F_BOA_rad)
        print("    F_rad^loss  = %.2f W m-2  " % F_loss)

        # Store previous atmosphere for reference
        atm_hist.append(copy.deepcopy(atm))

        # Plot
        plt_title = "Step %d:     $|dT_{comp}|$ = %.1f K     $|F_{rad}^{loss}|$ = %.1f W m$^{-2}$" % (step, dtmp_comp, F_loss)
        plot_radiative_eqm(atm_hist, atm_orig, plt_title)

        # Convergence (minimal temperature change between iterations => energy balance)
        success = flag_previous and (dtmp_comp < dtmp_conv) and (drel_dt < drel_dt_conv)
        flag_previous = (dtmp_comp < dtmp_conv)
        if flag_previous:
            print("    almost converged")

        step += 1
        print(" ")


    if not success:
        print("Stopping radiative equilibrium iterations without success")
    else:
        print("Found global radiative equilibrium")

    return atm

