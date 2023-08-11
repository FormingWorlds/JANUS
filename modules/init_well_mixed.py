import numpy as np
from utils.GeneralAdiabat import atm_z, molar_mass, cpv, p_sat, interpolate_atm
from utils.atmosphere_column import atmos
from utils.ClimateUtilities import integrator

def isothermal(lnP, lnT, atm):
    return 0

def no_condensation( atm, idx, wet_list, dry_list, prs_reset):

    if idx==0:
        
        # Renormalize mixing ratios to ensure sum == 1
        if sum(atm.vol_list.values()) != 1.:
            vol_list_new = {}
            for vol in atm.vol_list.keys():
                vol_list_new[vol]   = atm.vol_list[vol] / sum(atm.vol_list.values())
            for vol in atm.vol_list.keys():
                atm.vol_list[vol] = vol_list_new[vol]
        # Partial pressures scaled 
        for vol in atm.vol_list:
            atm.p_vol[vol][idx] = atm.vol_list[vol] * atm.p[idx]
            #atm.p_vol[vol][idx] = atm.vol_list[vol] * p_tot_pre_condensation
            # Mean gas phase molar mass
            atm.mu[idx]       += atm.vol_list[vol] * molar_mass[vol]
    
    else:  
        
        # Total pressure & partial pressures of condensing phases
        p_cond_sum = 0
        for vol in atm.vol_list:
            if vol in wet_list:
                atm.p_vol[vol][idx] = p_sat(vol, atm.tmp[idx])
                p_cond_sum += atm.p_vol[vol][idx]
                
        # Calculate the total partial pressure of dry species
        p_dry_tot = atm.p[idx] - p_cond_sum
        dry_frac_sum = 0
        # calculate sum of fractions of dry species for partial pressure calculation
        for vol in atm.vol_list:
            if vol in dry_list:
                dry_frac_sum += atm.x_gas[vol][idx-1]
        # Calculate the individual partial pressures of dry species
        for vol in atm.vol_list:
            if vol in dry_list:
                atm.p_vol[vol][idx] = p_dry_tot * (atm.x_gas[vol][idx-1]/dry_frac_sum)
                
       
    # Calculate mean molar mass
    atm.mu[idx]   = 0.
    
    for vol in atm.vol_list.keys():
        atm.mu[idx]   += molar_mass[vol] * atm.p_vol[vol][idx]
    atm.mu[idx] /= atm.p[idx]
    
    # Update condensate mixing ratios
    for vol in atm.vol_list.keys():
        atm.x_cond[vol][idx] = 0.
        atm.xc[idx]     += atm.x_cond[vol][idx]


    for vol in atm.vol_list.keys():    
        if idx == 0:
            atm.x_gas[vol][idx] = atm.vol_list[vol]
        else:
            atm.x_gas[vol][idx] =  ( 1-atm.xc[idx] ) * atm.p_vol[vol][idx] /  atm.p[idx]
        
        if vol in dry_list:
            atm.xd[idx]          += atm.x_gas[vol][idx]
        if vol in wet_list:
            atm.xv[idx]          += atm.x_gas[vol][idx]
        
        # Mean cp of both gas phase and retained condensates
        atm.cp[idx]    += atm.x_gas[vol][idx] * cpv(vol, atm.tmp[idx]) 

    # Renormalize cp w/ molar concentration (= cp_hat, Eq. 1, Li+2018)
    atm.cp[idx]  = atm.cp[idx] / ( atm.xd[idx] + atm.xv[idx] )

    # Dry concentration floor
    atm.xd[idx]  = np.amax([atm.xd[idx], 1e-10])
    
    # Loop through species to determine wet_list and dry_list for next level
    for vol in atm.vol_list.keys():
        
        if atm.vol_list[vol] > 0. and vol not in dry_list and vol not in wet_list:
            dry_list.append(vol)
    
    return atm, wet_list, dry_list


def ini_wm_iso(atm):

    new_psurf = 0
    new_p_vol = {}
    wet_list = []
    dry_list = []
    for vol in atm.vol_list.keys():
        if atm.vol_list[vol] * atm.ps > p_sat(vol, atm.ts):
            new_psurf += p_sat(vol,atm.ts)
            new_p_vol[vol] = p_sat(vol,atm.ts)
        else:
            new_psurf += atm.vol_list[vol] * atm.ps
            new_p_vol[vol] = atm.vol_list[vol] * atm.ps
            
    if new_psurf != atm.ps:
        Tsurf = atm.ts
        alpha = atm.alpha_cloud
        toa_heating = atm.toa_heating
        minT = atm.minT

        for vol in atm.vol_list.keys():
            atm.vol_list[vol] = new_p_vol[vol] / new_psurf

        atm = atmos(Tsurf, new_psurf, atm.ptop, atm.planet_radius, atm.planet_mass, vol_mixing=atm.vol_list, trppT=atm.trppT, minT=minT)
        
        atm.alpha_cloud = alpha
        atm.toa_heating = toa_heating
        
    for vol in atm.vol_list.keys():
        dry_list.append(vol)
    
    ### Initialization
    
    # Initialize the tuple solution
    moist_tuple     = [] #[tuple([np.log(atm.ps), atm.ts])] 
   
    # Negative increment to go from ps to ptop < ps       
    step            = -.01

    # Integration counter
    idx             = 0  
    
    # Surface condensation
    atm, wet_list, dry_list  = no_condensation(atm, idx, wet_list, dry_list, prs_reset=False)

    # Create the integrator instance                                              
    int_slope       = integrator(isothermal, np.log(atm.ps), np.log(atm.ts), step)

    # Update parameters used in the slope function dlntT/dlnP
    int_slope.setParams(atm)
    
    ### Integration of full general adiabat
    while atm.p[idx] > atm.ptop:
        
        # Execute the Runge-Kutta integrator, fill array of tuples
        moist_tuple.append(int_slope.next())

        # Fill next T,P values
        atm.p[idx+1]    = np.exp(int_slope.x)
        atm.tmp[idx+1]  = np.exp(int_slope.y)

        # Calculate next local gravity and height
        atm = atm_z(atm, idx)

        # print("RK4 step, idx:", idx, round(atm.p[idx+1],5), round(atm.tmp[idx+1],5))

        # Set next level to calculate
        idx             += 1
        atm.ifatm[idx]  = idx

        # Calculate condensation at next level
        atm, wet_list, dry_list    = no_condensation(atm, idx, wet_list, dry_list, prs_reset=False)

    # Interpolate
    atm = interpolate_atm(atm)
    
    # Rescale mixing ratios for plotting purposes
    atm.xd[:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    atm.xv[:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    for vol in atm.vol_list.keys():
        if atm.vol_list[vol] > 1e-10:
            atm.x_cond[vol][:] *= atm.alpha_cloud / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
            atm.x_gas[vol][:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    atm.xc[:] *= atm.alpha_cloud / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    
    return atm

