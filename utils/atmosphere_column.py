# atmosphere_column.py
# class for atmospheric column data

import numpy as np
from utils import phys

class atmos:
    """Atmosphere class    
    
    Stores compositional and thermodynamic information for the column.    
    Also stores planetary parameters.   

    """
    
    def __init__(self, T_surf: float, P_surf: float, P_top: float, pl_radius: float, pl_mass: float, 
                 vol_mixing: dict = {}, vol_partial: dict = {}, 
                 calc_cf: bool=False, 
                 trppT: float = 290):
        
        """Inits atmos class.

        One of either vol_mixing or vol_partial must be passed in. 
        If vol_partial is passed, then the value of P_surf is recalculated as the sum of partial pressures.

        Parameters
        ----------
            T_surf : float
                Surface temperature [K]
            P_surf : float
                Surface pressure [Pa]
            P_top : float
                Pressure at top of column [Pa]
            pl_radius : float
                Radius of rocky part of planet [m]
            pl_mass : float
                Mass of rocky part of planet [kg]

            vol_mixing : dict
                Dictionary of volatiles (keys) and mixing ratios (values)
            vol_partial: dict
                Dictionary of volatiles (keys) and partial pressures (values)

            calc_cf : bool
                Calculate contribution function?

            trppT : float
                Tropopause temperature
                
        """

        # Parse volatiles 
        if (len(vol_mixing) == 0) and (len(vol_partial) == 0):
            raise Exception("Either vol_mixing OR vol_partial must be passed to atmos.__init__ function!\nNeither were.")
        if (len(vol_mixing) > 0) and (len(vol_partial) > 0):
            raise Exception("Either vol_mixing OR vol_partial must be passed to atmos.__init__ function!\nBoth were.")
        
        if len(vol_mixing) > 0:  # Set by mixing ratio
            self.ps = P_surf
            if (P_surf <= 0.0):
                raise Exception("Surface pressure passed to atmos.__init__ must be positive!\nValue passed = '%g'" % P_surf)
            
            tot_mixing =  float(sum(vol_mixing.values()))  # Ensure mixing ratios add up to unity
            self.vol_list = {}
            for key in vol_mixing.keys():
                self.vol_list[key] = vol_mixing[key]/tot_mixing

        
        if len(vol_partial) > 0: # Set by partial pressure
            self.ps = float(sum(vol_partial.values()))
            self.vol_list = {}
            for key in vol_partial.keys():
                self.vol_list[key] = vol_partial[key]/self.ps

        # Initialise other variables
        self.alpha_cloud 	= 0.0 	    	# The fraction of condensate retained in the column; 1 -> Li et al 2018; 0 -> full rainout
        
        self.ts 			= T_surf		# Surface temperature, K

        self.ptop 			= P_top 			# Top pressure in Pa
        self.nlev 			= 10000  	   	# Number of vertical levels for adiabat integration
        self.step    		= 0.01  		# Adjust to match self.nlev
        self.nlev_save		= 100   		# Number of levels to save object
        self.p 				= np.zeros(self.nlev) 	   		# np.ones(self.nlev)
        self.pl 			= np.zeros(self.nlev+1)    		# np.ones(self.nlev+1)

        self.trppT          = trppT                 # Fixed value [K]
        self.trppidx		= 0 				 	# Tropopause: idx
        self.trppP 			= 0 				 	# Tropopause: prs

        self.dt 			= 0.5 							# days

        self.toa_heating    = 0. 							# W/m^2
        self.star_lum       = 0.0							# L_sun

        self.albedo_s   	= 0.0 							# surface albedo
        self.albedo_pl   	= 0.175 						# Bond albedo (scattering)
        self.zenith_angle  	= 54.55							# solar zenith angle, Hamano+15 (arccos(1/sqrt(3) = 54.74), Wordsworth+ 10: 48.19 (arccos(2/3)), see Cronin 14 (mu = 0.58 -> theta = arccos(0.58) = 54.55) for definitions

        self.planet_mass    = pl_mass
        self.planet_radius  = pl_radius
        self.grav_s 		= phys.G*self.planet_mass/(self.planet_radius**2) # m s-2
        
        self.tmp 			= np.zeros(self.nlev)      		# self.ts*np.ones(self.nlev)
        self.tmpl 			= np.zeros(self.nlev+1)
        self.Rcp    		= 2./7. 						# standard earth air
        self.n_species 		= 7
        self.mixing_ratios 	= np.zeros([self.n_species,self.nlev])
        
        # self.bands 			= np.concatenate((np.arange(0,3000,20),np.arange(3000,9000,50),np.arange(9000,24500,500))) # cm
        self.bands 			= np.concatenate((np.arange(0,3000,25),np.arange(3000,11000,50),np.arange(11000,30500,500))) # cm, 318 bands: HITEMP-compatible spacing
        
        self.band_centres 	= (self.bands[1:] + self.bands[:-1]) / 2
        self.band_widths 	= np.diff(self.bands)
        self.nbands 	    = np.size(self.bands)-1

        # Level-dependent quantities
        self.p_vol 			= {} # Gas phase partial pressures
        self.pl_vol 		= {} # Gas phase partial pressures
        self.x_gas 			= {} # Gas phase molar concentration
        self.x_cond         = {} # Condensed phase molar concentration
        self.grav_z			= np.zeros(self.nlev) # Local gravity
        self.z 				= np.zeros(self.nlev) # Atmospheric height
        self.mu 			= np.zeros(self.nlev) # Mean molar mass in level
        self.xd 			= np.zeros(self.nlev) # Molar concentration of dry gas
        self.xv 			= np.zeros(self.nlev) # Molar concentration of moist gas
        self.xc 			= np.zeros(self.nlev) # Molar concentration of condensed phase
        self.rho            = np.zeros(self.nlev) # Density of atmosphere at a given level
        self.ifatm 			= np.zeros(self.nlev) # Defines nth level to which atmosphere is calculated
        self.cp      		= np.zeros(self.nlev) # Mean heat capacity

        # Define T and P arrays from surface up
        self.tmp[0]         = self.ts         		# K
        self.p[0]           = self.ps         		# Pa
        self.z[0]           = 0         			# m
        self.grav_z[0]      = self.grav_s 			# m s-2


        # H2O floor to prevent NaNs
        self.vol_list["H2O"] = np.max( [ self.vol_list["H2O"], 1e-30 ] )

        # Instantiate object dicts and arrays
        for vol in self.vol_list.keys():
            # Instantiate as zero
            self.p_vol[vol]      = np.zeros(self.nlev)
            self.pl_vol[vol]     = np.zeros(self.nlev+1)
            self.x_gas[vol]      = np.zeros(self.nlev)
            self.x_cond[vol]     = np.zeros(self.nlev)
            
            #self.x_ocean[vol]    = 0.

            # Surface partial pressures
            self.p_vol[vol][0]   = self.ps * self.vol_list[vol]

        # Radiation heating and fluxes
        self.LW_flux_up 			= np.zeros(self.nlev)				# W/m^2
        self.LW_flux_down 			= np.zeros(self.nlev)				# W/m^2
        self.LW_flux_net			= np.zeros(self.nlev)				# W/m^2
        self.LW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev])	# W/m^2/(band)
        self.LW_heating				= np.zeros(self.nlev)				# K/day
        self.SW_flux_up 			= np.zeros(self.nlev)				# W/m^2
        self.SW_flux_down 			= np.zeros(self.nlev)				# W/m^2
        self.SW_flux_net			= np.zeros(self.nlev)				# W/m^2
        self.SW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev])	# W/m^2/(band)
        self.SW_heating				= np.zeros(self.nlev)				# K/day
        self.flux_up_total			= np.zeros(self.nlev)				# W/m^2
        self.flux_down_total		= np.zeros(self.nlev)				# W/m^2
        self.net_flux				= np.zeros(self.nlev)				# W/m^2
        self.net_spectral_flux	 	= np.zeros([self.nbands,self.nlev])	# W/m^2/(band)
        self.net_heating 			= np.zeros(self.nlev) 				# K/day
        # Contribution function arrays
        if calc_cf == True:
            self.cff 					= np.zeros(self.nlev) 				# normalised
            self.cff_i					= np.zeros([self.nbands,self.nlev]) # cf per band
            self.LW_flux_up_i 			= np.zeros([self.nbands,self.nlev])

    def write_PT(self,filename: str="output/PT.tsv", punit:str = "Pa"):
        """Write PT profile to file, with descending pressure.

        Useful for when run outside of PROTEUS and for debugging.
        Designed with VULCAN compatibility in mind.

        Parameters
        ----------
            filename : string
                Output filename
            punit : string
                Pressure unit to use. Options: Pa, bar, dyne/cm2, atm.

        """

        p_scalefactor = 1.0
        match punit:
            case "Pa":
                p_scalefactor = 1.0
            case "bar":
                p_scalefactor = 1.0e-5
            case "dyne/cm2":
                p_scalefactor = 1.0e5
            case "atm":
                p_scalefactor = 1.01325e5
            case _:
                raise Exception("Unrecognised pressure unit '%s'"%punit)

        p_save = np.array(self.p) * p_scalefactor
        T_save = np.array(self.tmp)

        X = np.array([p_save,T_save]).T[::-1]

        header = '# (%s)\t(K)\nPressure\tTemp' % punit

        np.savetxt(filename,X,fmt='%1.5e',header=header,comments='',delimiter='\t')


