# atmosphere_column.py
# class for atmospheric column data

import numpy as np

surface_pressure 	= 1e+5 						# Pa
top_pressure 		= 1e-7*surface_pressure 	# Pa
n_vertical_levels 	= 10000  					# For computation
timestep 			= 0.5
n_absorbing_species = 7
n_bands 			= 300

R_universal = 8.31446261815324 # Universal gas constant, J.K-1.mol-1

class atmos:
	'''
	Atmosphere class
	'''
	def __init__(self):
		self.ps 			= surface_pressure 	 	   # Surface pressure in Pa
		self.ptop 			= top_pressure 		 	   # Top pressure in Pa
		self.nlev 			= n_vertical_levels  	   # Number of vertical levels
		self.p 				= np.zeros(self.nlev) 	   # np.ones(self.nlev)
		self.pl 			= np.zeros(self.nlev-1)    # np.ones(self.nlev+1)
		self.dt 			= timestep
		self.ts 			= 300.0
		self.tmp 			= np.zeros(self.nlev)      # self.ts*np.ones(self.nlev)
		self.tmpl 			= np.zeros(self.nlev-1)
		self.Rcp 			= 2./7.
		self.n_species 		= n_absorbing_species
		self.mixing_ratios 	= np.zeros([self.n_species,self.nlev])
		self.fluxes 		= self.atmos_fluxes(self.nlev)
		self.bands 			= np.concatenate((np.arange(0,3000,20),np.arange(3000,9000,50),np.arange(9000,24500,500)))
		self.band_centres 	= (self.bands[1:] + self.bands[:-1]) / 2
		self.band_widths 	= np.diff(self.bands)

		# Species-dependent quantities
		self.p_vol 			= {} # Gas phase partial pressures
		self.x_gas 			= {} # Gas phase molar concentration
		self.x_gasl 		= {} # Gas phase molar concentration staggered
		self.x_cond         = {} # Condensed phase molar concentration
		self.mr_gas 		= {} # Gas phase molar mixing ratio (relative to all gas)
		self.mr_cond        = {} # Condensed phase molar mixing ratio (relative to all gas)
		
		# Level-dependent quantities
		self.xd 			= np.zeros(self.nlev)	# Molar concentration of dry gas
		self.xv 			= np.zeros(self.nlev)	# Molar concentration of moist gas
		self.xc 			= np.zeros(self.nlev)	# Molar concentration of condensed phase
		self.mrd 			= np.zeros(self.nlev)	# Molar mixing ratio of 'dry' gas (relative to gas phase)
		self.mrv 			= np.zeros(self.nlev)	# Molar mixing ratio of 'condensing' gas (relative to gas)
		self.mrc			= np.zeros(self.nlev)	# Molar mixing ratio of cloud phase (relative to gas)
		self.vol_list 		= []			     	# names of all species present
		self.ifatm 			= np.zeros(self.nlev) 	# Defines n level to which atmosphere is calculated
		self.cp      		= np.zeros(self.nlev)   # Heat capactiy depending on molar concentration ratio
		self.cp_mr     		= np.zeros(self.nlev)   # Heat capactiy depending on mixing ratio

	class atmos_fluxes:
		'''
		Fluxes class
		'''
		def __init__(self,nlev):
			self.nlev_flux 				= nlev
			self.LW_flux_up 			= np.zeros(nlev)
			self.LW_spectral_flux_up 	= np.zeros([n_bands,nlev])
			self.total_heating 			= np.zeros(nlev) 				# K s^-1
