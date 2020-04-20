# atmosphere_column.py
# class for atmospheric column data

import numpy as np

R_universal = 8.31446261815324 # Universal gas constant, J.K-1.mol-1

class atmos:
	'''
	Atmosphere class
	'''
	def __init__(self, T_surf, P_surf, vol_list):
		self.ps 			= P_surf 	 	# Surface pressure, Pa
		self.ts 			= T_surf		# Surface temperature, K
		self.vol_list 		= vol_list		# Names + mixing ratios dict

		self.ptop 			= 1 			# Top pressure in Pa
		self.nlev 			= 10000  	   	# Number of vertical levels for adiabat integration
		self.nlev_save		= 100   		# Number of levels to save object
		self.p 				= np.zeros(self.nlev) 	   		# np.ones(self.nlev)
		self.pl 			= np.zeros(self.nlev+1)    		# np.ones(self.nlev+1)

		self.trpp 			= np.zeros(3) 				 	# Tropopause: idx, prs, tmp
		
		self.dt 			= 0.5 							# days

		self.toa_heating    = 0. 							# W/m^2
		self.star_lum       = 0. 							# L_sun

		self.albedo_s   	= 0.1 							# surface albedo
		self.albedo_pl   	= 0.2 							# planetary albedo (scattering)
		self.zenith_angle  	= 38							# solar zenith angle
		
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

		# Species-dependent quantities
		self.p_vol 			= {} # Gas phase partial pressures
		self.x_gas 			= {} # Gas phase molar concentration
		self.x_cond         = {} # Condensed phase molar concentration
		self.mr_gas 		= {} # Gas phase molar mixing ratio (relative to all gas)
		self.mr_cond        = {} # Condensed phase molar mixing ratio (relative to all gas)
		self.x_ocean		= {} # Surface condensed 'overpressure'
		
		# Level-dependent quantities
		self.xd 			= np.zeros(self.nlev)	# Molar concentration of dry gas
		self.xv 			= np.zeros(self.nlev)	# Molar concentration of moist gas
		self.xc 			= np.zeros(self.nlev)	# Molar concentration of condensed phase
		self.mrd 			= np.zeros(self.nlev)	# Molar mixing ratio of 'dry' gas (relative to gas phase)
		self.mrv 			= np.zeros(self.nlev)	# Molar mixing ratio of 'condensing' gas (relative to gas)
		self.mrc			= np.zeros(self.nlev)	# Molar mixing ratio of cloud phase (relative to gas)
		self.ifatm 			= np.zeros(self.nlev) 	# Defines n level to which atmosphere is calculated
		self.cp      		= np.zeros(self.nlev)   # Heat capacity depending on molar concentration ratio
		self.cp_mr     		= np.zeros(self.nlev)   # Heat capacity depending on mixing ratio

		# Define T and P arrays from surface up
		self.tmp[0]         = T_surf         		# K
		self.p[0]           = P_surf         		# Pa


		# H2O floor to prevent NaNs
		self.vol_list["H2O"] = np.max( [ self.vol_list["H2O"], 1e-30 ] )

		# Instantiate object dicts and arrays
		for vol in self.vol_list.keys():
		    # Instantiate as zero
		    self.p_vol[vol]      = np.zeros(self.nlev)
		    self.x_gas[vol]      = np.zeros(self.nlev)
		    self.x_cond[vol]     = np.zeros(self.nlev)
		    self.mr_gas[vol]     = np.zeros(self.nlev)
		    self.mr_cond[vol]    = np.zeros(self.nlev)
		    self.x_ocean[vol]    = 0.

		    # Surface partial pressures
		    self.p_vol[vol][0]   = self.ps * vol_list[vol]

		    

		# Radiation heating and fluxes
		self.LW_flux_up 			= np.zeros(self.nlev)				# W/m^2
		self.LW_flux_down 			= np.zeros(self.nlev)				# W/m^2
		self.LW_flux_net			= np.zeros(self.nlev)				# W/m^2
		self.LW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev])		# W/m^2/(band)
		self.LW_heating				= np.zeros(self.nlev)				# K/day
		self.SW_flux_up 			= np.zeros(self.nlev)				# W/m^2
		self.SW_flux_down 			= np.zeros(self.nlev)				# W/m^2
		self.SW_flux_net			= np.zeros(self.nlev)				# W/m^2
		self.SW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev])		# W/m^2/(band)
		self.SW_heating				= np.zeros(self.nlev)				# K/day
		self.flux_up_total			= np.zeros(self.nlev)				# W/m^2
		self.flux_down_total		= np.zeros(self.nlev)				# W/m^2
		self.net_flux				= np.zeros(self.nlev)				# W/m^2
		self.net_spectral_flux	 	= np.zeros([self.nbands,self.nlev])		# W/m^2/(band)
		self.net_heating 			= np.zeros(self.nlev) 				# K/day
		