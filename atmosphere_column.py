# atmosphere_column.py
# class for atmospheric column data
# MDH 07/05/19

import numpy as np

surface_pressure = 1.e3 # Surface pressure in mbar
top_pressure = 5.0
n_vertical_levels = 50
timestep = 0.5
n_absorbing_species = 2
n_bands = 9

class atmos:
	'''
	Atmosphere class
	'''
	def __init__(self):
		self.ps = surface_pressure # Surface pressure in mbar
		self.ptop = top_pressure # Top pressure in mbar
		self.nlev = n_vertical_levels # Number of vertical levels
		self.p = np.ones(self.nlev)
		self.pl = np.ones(self.nlev+1)
		self.dt = timestep
		self.temp = 300.0*np.ones(self.nlev)
		self.templ = 300.0*np.ones(self.nlev+1)
		self.ts = 300.0
		self.Rcp = 2./7.
		self.n_species = n_absorbing_species
		self.mixing_ratios = np.zeros([self.n_species,self.nlev])
		self.fluxes = self.atmos_fluxes(self.nlev)
# TODO: read in bands from spectral file
		self.bands = 0.01/np.array([1.000000000E-02,2.500000000E-05,1.818181818E-05,1.250000000E-05,1.333333333E-05,8.333333333E-06,8.928571429E-06,7.518796992E-06,6.666666667E-06,3.338898164E-06])
		self.band_centres = (self.bands[1:] + self.bands[:-1]) / 2
		self.band_widths = np.diff(self.bands)

	class atmos_fluxes:
		'''
		Fluxes class
		'''
		def __init__(self,nlev):
			self.nlev_flux = nlev
			self.LW_flux_up = np.zeros(nlev)
			self.LW_spectral_flux_up = np.zeros([n_bands,nlev])
			self.total_heating = np.zeros(nlev)
