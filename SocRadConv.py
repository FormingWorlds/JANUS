'''
MDH 28/01/19

Socrates radiative-convective model
'''

import numpy as np
import math,phys
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib
import SocRadModel
from atmosphere_column import atmos


# # Font settings
# import matplotlib.pylab as pylab
# import matplotlib.font_manager as fm
# font = fm.FontProperties(family = 'Helvetica', fname = '/Users/tim/Dropbox/work/matplotlib_fonts/Helvetica/Helvetica.ttf')

"""------------------------------------------------------------------------ """
"""------------------------Thermodynamic constants------------------------- """
"""------------------------------------------------------------------------ """

R = phys.water.R                         # J/(kg*K) specific gas constant of water vapor
Rcp = phys.water.Rcp                     # cp in J/(kg*K) specific heat constant of water vapor
L = phys.water.L_vaporization            # J/kg, latent heat of condensation of water vapor at 300K
esat = phys.satvps_function(phys.water)  # Saturation vapor pressure, arguments=T,T0,e0,MolecularWeight,LatentHeat
Tref = 350.                              # Reference temperature
pref = esat(Tref)                        # Reference pressure

#Dew point temperature
def Tdew(p):
    return Tref/(1-(Tref*R/L)*math.log(p/pref))

""" Moist adjustment switch """
Moist_Adjustment = True

def surf_Planck_nu(atm):
    h = 6.63e-34
    c = 3.0e8
    kb = 1.38e-23
    B = np.zeros(len(atm.band_centres))
    c1 = 1.191042e-5
    c2 = 1.4387752
    for i in range(len(atm.band_centres)):
        nu = atm.band_centres[i]
        B[i] = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))

    B = B * atm.band_widths/1000.0
    return B

def RadConvEqm(output_dir, time_current, Tg, stellar_toa_heating, p_s, h2o_ratio, co2_ratio, h2_ratio, ch4_ratio, co_ratio, n2_ratio, o2_ratio, he_ratio):
    #--------------------Set radmodel options-------------------
    #---Instantiate the radiation model---

    atm = atmos()

    #---Set up pressure array (a global)----
    atm.ps = p_s
    pstart = .995*atm.ps
    rat = (atm.ptop/pstart)**(1./atm.nlev)
    logLevels = [pstart*rat**i for i in range(atm.nlev+1)]
    logLevels.reverse()
    levels = [atm.ptop + i*(pstart-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl = np.array(logLevels)
    atm.p = (atm.pl[1:] + atm.pl[:-1]) / 2


    #==============Now do the calculation====================================

    atm.ts = Tg
    atm.Rcp = 2./7.
    atm.temp = atm.ts*(atm.p/atm.p[-1])**atm.Rcp  #Initialize on an adiabat
    atm.temp  = np.where(atm.temp<atm.ts/1.5,atm.ts/1.5,atm.temp)
    # atm.n_species = 2
    atm.n_species = 7
    Moist_adiabat=[Tdew(pp) for pp in atm.p]
    
    # # Water vapour
    # atm.mixing_ratios[0] = 1.e-5
    # # CO2
    # atm.mixing_ratios[1] = 1.e-5

    atm.mixing_ratios[0] = h2o_ratio # H2O
    atm.mixing_ratios[1] = co2_ratio # CO2
    atm.mixing_ratios[2] = h2_ratio  # H2
    atm.mixing_ratios[3] = ch4_ratio # CH4
    atm.mixing_ratios[4] = co_ratio  # CO
    atm.mixing_ratios[5] = n2_ratio  # N2
    atm.mixing_ratios[6] = o2_ratio  # O2


    # Initialise previous OLR and TOA heating to zero
    PrevOLR = 0.
    PrevMaxHeat = 0.
    PrevTemp = 0.*atm.temp[:]

    #---------------------------------------------------------
    #--------------Initializations Done-----------------------
    #--------------Now do the time stepping-------------------
    #---------------------------------------------------------
    matplotlib.rc('axes',edgecolor='k')
    for i in range(0,1):

        atm = steps(atm, stellar_toa_heating)

        #hack!
        # atm.temp[0] = atm.temp[1]

        if i % 5 == 0:
            if 1==2:
                plt.figure(figsize=(7,4))
                plt.semilogy(atm.temp,atm.p)
                plt.gca().invert_yaxis()
                plt.ylabel('Pressure [mb]')
                plt.xlabel('Temperature [K]')
                plt.gca().xaxis.label.set_color('white')
                plt.tick_params(axis='x', colors='white')
                plt.gca().yaxis.label.set_color('white')
                plt.tick_params(axis='y', colors='white')
                # plt.show()
            #print("OLR " + str(atm.LW_flux_up[-1]))
            #print("OLR change " + str(atm.LW_flux_up[-1]-PrevOLR))
            # print("Max heating " + str(np.max(atm.total_heating)))
            #print("Max dT " + str(abs(np.max(atm.temp-PrevTemp[:]))))
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,7))
            ax1.semilogy(atm.temp,atm.p)
            ax1.semilogy(Moist_adiabat,atm.p,'r',label='Moist adiabat')
            ax1.invert_yaxis()
            ax1.set_xlabel('Temperature [K]')
            ax1.set_ylabel('Pressure [mb]')
            ax2.plot(atm.band_centres,atm.LW_spectral_flux_up[:,0]/atm.band_widths,'k')
#            ax2.plot(atm.band_centres,surf_Planck_nu(atm)/atm.band_widths,'k--')
            ax2.set_xlim([0,8*atm.ts])
            ax2.set_ylabel('Spectral Flux')
            ax2.set_xlabel('Wavenumber')
            ax2.set_title('Spectral OLR')
            # plt.show()
            plt.savefig(output_dir+'/TP_profile_'+str(round(time_current))+'.png', bbox_inches="tight")
            plt.show()
            plt.close(fig)
            print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)))

        # Reduce timestep if heating not converging
        if abs(np.max(atm.temp-PrevTemp[:])) < 0.05 or abs(atm.temp[0]-atm.temp[1]) > 3.0:
            print("reducing timestep")
            atm.dt  = atm.dt*0.99


        #if abs(atm.LW_flux_up[-1]-PrevOLR) < 10. or abs(np.max(atm.temp-PrevTemp[:])) < 10.:
        #print(str(atm.LW_flux_up[0])+","+str(PrevOLR)+","+str(abs(atm.LW_flux_up[-1]-PrevOLR)))
        if abs(atm.LW_flux_up[0]-PrevOLR) < 50. and i > 5 :
           print("break")
           #print(PrevTemp[:]-atm.temp)
           break    # break here

        PrevOLR = atm.LW_flux_up[0]
        PrevMaxHeat = abs(np.max(atm.total_heating))
        PrevTemp[:] = atm.temp[:]

    # # plot equilibrium temperature profile
    # plt.figure()
    # plt.semilogy(atm.temp,atm.p)
    # plt.gca().invert_yaxis()
    # plt.ylabel('Pressure [mb]')
    # plt.xlabel('Temperature [K]')
    # plt.savefig(output_dir+'/T_profile_'+str(round(time_current))+'.pdf', bbox_inches="tight")

    # Write TP and spectral flux profiles for later plotting
    out_a = np.column_stack( ( atm.temp, atm.p ) )
    np.savetxt( output_dir+str(int(time_current))+"_atm_TP_profile.dat", out_a )
    out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths ) )
    np.savetxt( output_dir+str(int(time_current))+"_atm_spectral_flux.dat", out_a )

    return atm.LW_flux_up[-1], atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths

#--------- Importing thermodynamical properties of gases -----------

def L_heat(switch): # Molar latent heat for gases considered, in J.mol-1
    if switch == 'H2O':
        return phys.water.L_vaporization*phys.water.MolecularWeight*1e-3 # Conversion from J.kg-1 to J.mol-1 (molecular weight is in g/mol in phys.f90)
    if switch == 'CH4':
        return phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
    if switch == 'CO2':
        return phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
    if switch == 'CO':
        return phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
    if switch == 'N2':
        return phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
    if switch == 'O2':
        return phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
    if switch == 'H2':
        return phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
    if switch == 'He':
        return phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
    if switch == 'NH3':
        return phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3

R_universal=8.31446261815324 # Universal gas constant, J.K-1.mol-1, should probably go elsewhere since it's a constant

def cpv(switch): # Molar heat capacities for gases considered, in J.K-1.mol-1
    if switch == 'H2O':
        return phys.water.cp*phys.water.MolecularWeight*1e-3
    if switch == 'CH4':
        return phys.methane.cp*phys.methane.MolecularWeight*1e-3
    if switch == 'CO2':
        return phys.co2.cp*phys.co2.MolecularWeight*1e-3
    if switch == 'CO':
        return phys.co.cp*phys.co.MolecularWeight*1e-3
    if switch == 'N2':
        return phys.n2.cp*phys.n2.MolecularWeight*1e-3
    if switch == 'O2':
        return phys.o2.cp*phys.o2.MolecularWeight*1e-3
    if switch == 'H2':
        return phys.h2.cp*phys.h2.MolecularWeight*1e-3
    if switch == 'He':
        return phys.he.cp*phys.he.MolecularWeight*1e-3
    if switch == 'NH3':
        return phys.nh3.cp*phys.nh3.MolecularWeight*1e-3
    
def cpc(switch): # Molar heat capacities for condensed phases considered. J.K-1.mol-1. # ToDo: replace with liquid or solid phase values
    if switch == 'H2O':
        return phys.water.cp*phys.water.MolecularWeight*1e-3 
    if switch == 'CH4':
        return phys.methane.cp*phys.methane.MolecularWeight*1e-3
    if switch == 'CO2':
        return phys.co2.cp*phys.co2.MolecularWeight*1e-3
    if switch == 'CO':
        return phys.co.cp*phys.co.MolecularWeight*1e-3
    if switch == 'N2':
        return phys.n2.cp*phys.n2.MolecularWeight*1e-3
    if switch == 'O2':
        return phys.o2.cp*phys.o2.MolecularWeight*1e-3
    if switch == 'H2':
        return phys.h2.cp*phys.h2.MolecularWeight*1e-3
    if switch == 'He':
        return phys.he.cp*phys.he.MolecularWeight*1e-3
    if switch == 'NH3':
        return phys.nh3.cp*phys.nh3.MolecularWeight*1e-3

def xv(switch): # Molar abundances of vapor phase/mole of gas mixture 
    if switch == 'H2O':
        return 0.5
    if switch == 'CH4':
        return 0.
    if switch == 'CO2':
        return 0.
    if switch == 'CO':
        return 0.
    if switch == 'N2':
        return 0.
    if switch == 'O2':
        return 0.
    if switch == 'H2':
        return 0.
    if switch == 'He':
        return 0.
    if switch == 'NH3':
        return 0.
    
def xc(switch): # Molar abundances of condensed phase/mole of gas mixture
    if switch == 'H2O':
        return 0.
    if switch == 'CH4':
        return 0.
    if switch == 'CO2':
        return 0.
    if switch == 'CO':
        return 0.
    if switch == 'N2':
        return 0.
    if switch == 'O2':
        return 0.
    if switch == 'H2':
        return 0.
    if switch == 'He':
        return 0.
    if switch == 'NH3':
        return 0.
    
#function that returns dlnT/dlnP
def General_moist_adiabat(T,t):
    xd=0.5 # Molar abundance of non-condensible species/mole of gas mixture
    cpd=1000.*phys.air.MolecularWeight*1.e-3 # molar heat capacity of the dry species
    xv_cpv=xv('H2O')*cpv('H2O')+xv('CO2')*cpv('CO2') # add other species as needed (generalizable sum?)
    xc_cpc=xc('H2O')*cpc('H2O')+xc('CO2')*cpc('CO2')
    first_term=(xv('H2O')/xd)*(L_heat('H2O')/(R_universal*T))**2 + (xv('CO2')/xd)*(L_heat('CO2')/(R_universal*T))**2
    quadr_term=((xv('H2O')/xd)*L_heat('H2O')/(R_universal*T) + (xv('CO2')/xd)*L_heat('CO2')/(R_universal*T))**2
    sum_abundances=xv('H2O')+xv('CO2')
    sum_ratio_abundances=(xv('H2O')/xd)+(xv('CO2')/xd)

    num=(1+xv('H2O')*(L_heat('H2O')/(xd*R_universal*T)))
    denom=(1/R_universal)*(xd*cpd+xv_cpv+xc_cpc/(xd+sum_abundances)) + (first_term+quadr_term/(1.+sum_ratio_abundances))
    dTdt=num/denom
    return dTdt

# Dry adjustment routine
def dryAdj(atm):
    T = atm.temp
    p = atm.p
    #Rcp is a global
    #Downward pass
    for i in range(len(T)-1):
        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        pfact = (p1/p2)**atm.Rcp
        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) #Equal layer masses
                              #Not quite compatible with how
                              #heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            atm.temp[i] = T1
            atm.temp[i+1] = T2
    #Upward pass
    for i in range(len(T)-2,-1,-1):
        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        pfact = (p1/p2)**atm.Rcp
        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) #Equal layer masses
                              #Not quite compatible with how
                              #heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            atm.temp[i] = T1
            atm.temp[i+1] = T2            

#Moist adjustment routine.
def moistAdj(atm):

    T = atm.temp
    p = atm.p

    if Moist_Adjustment: # switch
        #initial condition
        T0=300.
        Tdew=odeint(General_moist_adiabat,T0,np.log(p)) # integrate dlnT/dlnP

        for i in range(len(T)-1): # Downward pass
            if T[i]<Tdew(p[i]):
                T[i]=Tdew(p[i]) # temperature stays the same during the phase change
        for i in range(len(T)-2,-1,-1): # Upward pass
            if T[i]<Tdew(p[i]):
                T[i]=Tdew(p[i])

#Define function to do time integration for n steps
def steps(atm, stellar_toa_heating):
    atm = SocRadModel.radCompSoc(atm, stellar_toa_heating)
    dT = atm.total_heating*atm.dt
    #Limit the temperature change per step
    dT = np.where(dT>5.,5.,dT)
    dT = np.where(dT<-5.,-5.,dT)
    #Midpoint method time stepping
    #changed call to r.  Also modified to hold Tg fixed
    atm = SocRadModel.radCompSoc(atm, stellar_toa_heating)
    dT = atm.total_heating*atm.dt
    #Limit the temperature change per step
    dT = np.where(dT>5.,5.,dT)
    dT = np.where(dT<-5.,-5.,dT)
    atm.temp += dT
    #
    dTmax = max(abs(dT)) #To keep track of convergence

    #   Do the surface balance
    kturb = .1
    atm.temp[-1] += -atm.dt*kturb*(atm.temp[-1] - atm.ts)
    #Dry adjustment step
    for iadj in range(10):
        dryAdj(atm)
        moistAdj(atm)
    Tad = atm.temp[-1]*(atm.p/atm.p[-1])**atm.Rcp
    #** Temporary kludge to keep stratosphere from getting too cold
    atm.temp = np.where(atm.temp<50.,50.,atm.temp)  #**KLUDGE
    #
    #Dummies for separate LW and stellar. **FIX THIS**
    fluxStellar = fluxLW = heatStellar = heatLW = np.zeros(atm.nlev)
    return atm
