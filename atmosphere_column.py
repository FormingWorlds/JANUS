# atmosphere_column.py
# class for atmospheric column data
# MDH 07/05/19

import numpy as np

surface_pressure 	= 1e5 						# Pa
top_pressure 		= 1e-7*surface_pressure 	# Pa
n_vertical_levels 	= 100
timestep 			= 0.5
n_absorbing_species = 7
n_bands 			= 300

class atmos:
	'''
	Atmosphere class
	'''
	def __init__(self):
		self.ps 			= surface_pressure 		# Surface pressure in Pa
		self.ptop 			= top_pressure 			# Top pressure in Pa
		self.nlev 			= n_vertical_levels 	# Number of vertical levels
		self.p 				= np.ones(self.nlev)
		self.pl 			= np.ones(self.nlev+1)
		self.dt 			= timestep
		self.ts 			= 300.0
		self.temp 			= self.ts*np.ones(self.nlev)
		self.templ 			= self.ts*np.ones(self.nlev+1)
		self.Rcp 			= 2./7.
		self.n_species 		= n_absorbing_species
		self.mixing_ratios 	= np.zeros([self.n_species,self.nlev])
		self.fluxes 		= self.atmos_fluxes(self.nlev)
		self.bands 			= np.concatenate((np.arange(0,3000,20),np.arange(3000,9000,50),np.arange(9000,24500,500)))
		self.band_centres 	= (self.bands[1:] + self.bands[:-1]) / 2
		self.band_widths 	= np.diff(self.bands)
		self.p_vol 			= {}
		self.vol_dry        = "H2"

	class atmos_fluxes:
		'''
		Fluxes class
		'''
		def __init__(self,nlev):
			self.nlev_flux 				= nlev
			self.LW_flux_up 			= np.zeros(nlev)
			self.LW_spectral_flux_up 	= np.zeros([n_bands,nlev])
			self.total_heating 			= np.zeros(nlev)
