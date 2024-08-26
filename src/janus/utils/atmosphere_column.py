# atmosphere_column.py
# class for atmospheric column data

import toml
import numpy as np
import netCDF4 as nc
from janus.utils import phys
from janus.utils.height import AtmosphericHeight
import os, copy, platform, shutil
import pwd

class atmos:
    
    def __init__(self, T_surf: float, P_surf: float, P_top: float, pl_radius: float, pl_mass: float, 
                 band_edges:list, vol_mixing: dict = {}, vol_partial: dict = {},
                 req_levels: int = 100, water_lookup: bool=False, alpha_cloud:float=0.0,
                 trppT: float = 290.0, minT: float = 1.0, maxT: float = 9000.0, do_cloud: bool=False,
                 re: float=0., lwm: float=0., clfr: float=0.,
                 albedo_s: float=0.0, albedo_pl: float=0.175, zenith_angle: float=54.74
                 ):
        
        """Atmosphere class    
    
        Stores compositional and thermodynamic information for the column.    
        Also stores planetary parameters.   

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
            band_edges : list
                List of band edges in nm, ascending

            vol_mixing : dict
                Dictionary of volatiles (keys) and mixing ratios (values)
            vol_partial: dict
                Dictionary of volatiles (keys) and partial pressures (values)

            req_levels : int
                Requested number of vertical levels
            alpha_cloud : float 
                Condensate retention fraction (1 -> Li et al 2018; 0 -> full rainout)
            water_lookup : bool
                Use lookup table for water thermodynamic values (e.g. L, c_p)
            trppT : float
                Tropopause temperature
            minT : float
                Temperature floor
            maxT : float
                Temperature ceiling
            re : float
                Effective radius of cloud droplets [m]
            lwm : float
                Liquid water mass fraction [kg/kg]
            clfr : float
                Water cloud fraction [adimensional]
            albedo_s : float
                Surface albedo
            albedo_pl : float
                Bond albedo (scattering) applied to self.toa_heating in socrates.py
            zenith_angle : float
                solar zenith angle, Hamano+15 (arccos(1/sqrt(3) = 54.74), Wordsworth+10: 48.19
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

        # Required volatiles
        required_vols = {"H2O","CO2","N2"}
        if len(required_vols.intersection(self.vol_list.keys())) < len(required_vols):
            raise Exception("Missing required volatiles!\nRequired vols = %s" % str(required_vols))
        
        # H2O floor to prevent NaNs
        self.vol_list["H2O"] = np.max( [ self.vol_list["H2O"], 1e-20 ] )
        
        # Remove volatiles with mixing ratio of zero
        # vol_list_old = copy.deepcopy(self.vol_list)
        # self.vol_list = {}
        # for key in vol_list_old.keys():
        #     if (vol_list_old[key] > 1.0e-30) or (key in required_vols):
        #         self.vol_list[key] = vol_list_old[key]


        # Initialise other variables
        self.alpha_cloud 	= alpha_cloud 	    	# The fraction of condensate retained in the column
        
        self.ts 			= T_surf		# Surface temperature, K

        self.ptop 			= P_top 			# Top pressure in Pa
        self.nlev 			= 10000  	   	# Number of vertical levels for adiabat integration
        self.step    		= 0.01  		# Adjust to match self.nlev
        self.nlev_save		= int(max(req_levels, 10))  		# Number of levels to save object
        self.p 				= np.zeros(self.nlev) 	   		# np.ones(self.nlev)
        self.pl 			= np.zeros(self.nlev+1)    		# np.ones(self.nlev+1)484854

        self.trppT          = trppT                 # Fixed value [K]
        self.trppidx		= 0 				 	# Tropopause: idx
        self.trppP 			= 0 				 	# Tropopause: prs
        self.minT           = minT                  # Minimum temperature allowed [K]
        self.maxT           = maxT                  # Maximum ^

        self.dt 			= 0.5 							# days

        self.instellation   = 0. 							# Instellation at planet's orbital separation, W/m^2
        self.toa_heating    = 0.                            # ASF

        self.inst_sf        = 3.0/8.0                       # Scale factor applied to instellation (see Cronin+14 for definitions)
        self.albedo_s   	= albedo_s                      # surface albedo
        self.albedo_pl   	= albedo_pl                     # Bond albedo (scattering) applied to self.toa_heating in socrates.py
        self.zenith_angle  	= zenith_angle                  # solar zenith angle, Hamano+15 (arccos(1/sqrt(3) = 54.74), Wordsworth+10: 48.19 (arccos(2/3)), see Cronin+14 for definitions

        self.planet_mass    = pl_mass
        self.planet_radius  = pl_radius
        self.grav_s 		= phys.G*self.planet_mass/(self.planet_radius**2) # m s-2
        
        self.tmp 			= np.zeros(self.nlev)      		# self.ts*np.ones(self.nlev)
        self.tmpl 			= np.zeros(self.nlev+1)
        self.Rcp    		= 2./7. 						# standard earth air
        self.n_species 		= 16
        self.mixing_ratios 	= np.zeros([self.n_species,self.nlev])

        self.water_lookup   = water_lookup
        
        self.tmp_magma      = 3000.0
        self.skin_d         = 0.01 # m
        self.skin_k         = 2.0  # W m-1 K-1
        
        # Level-dependent quantities
        self.p_vol 			= {} # Gas phase partial pressures
        self.pl_vol 		= {} # Gas phase partial pressures
        self.x_gas 			= {} # Gas phase molar concentration
        self.x_cond         = {} # Condensed phase molar concentration
        self.grav_z			= np.zeros(self.nlev) # Local gravity
        self.z 				= np.zeros(self.nlev)   # Atmospheric height (centres)
        self.zl 	    	= np.zeros(self.nlev+1) # Atmospheric height (edges)
        self.mu 			= np.zeros(self.nlev) # Mean molar mass in level
        self.xd 			= np.zeros(self.nlev) # Molar concentration of dry gas
        self.xv 			= np.zeros(self.nlev) # Molar concentration of moist gas
        self.xc 			= np.zeros(self.nlev) # Molar concentration of condensed phase
        self.rho            = np.zeros(self.nlev) # Density of atmosphere at a given level
        self.ifatm 			= np.zeros(self.nlev) # Defines nth level to which atmosphere is calculated
        self.cp      		= np.zeros(self.nlev) # Mean heat capacity
        self.height_error   = True # error when doing hydrostatic integration?

        # Define T and P arrays from surface up
        self.tmp[0]         = self.ts         		# K
        self.tmpl[0]         = self.ts         		# K
        self.p[0]            = self.ps         		# Pa
        self.pl[0]           = self.ps         		# Pa
        self.z[0]           = 0         			# m
        self.zl[0]
        self.grav_z[0]      = self.grav_s 			# m s-2


        # Instantiate object dicts and arrays
        for vol in self.vol_list.keys():
            # Instantiate as zero
            self.p_vol[vol]      = np.zeros(self.nlev)
            self.pl_vol[vol]     = np.zeros(self.nlev+1)
            self.x_gas[vol]      = np.zeros(self.nlev)
            self.x_cond[vol]     = np.zeros(self.nlev)
        
            # Surface partial pressures
            self.p_vol[vol][0]   = self.ps * self.vol_list[vol]


        # Spectral file
        self.has_contfunc = False
        self.band_edges = np.array(band_edges)  # units of [nm]
        self.bands_set = bool(len(band_edges) > 1)
        if self.bands_set:
            self.nbands 	    = np.size(self.band_edges)-1
            self.band_centres 	= (self.band_edges[1:] + self.band_edges[:-1]) / 2
            self.band_widths 	= np.abs(np.diff(self.band_edges))

            # Radiation heating and fluxes
            self.LW_flux_up 			= np.zeros(self.nlev_save)				# W/m^2
            self.LW_flux_down 			= np.zeros(self.nlev_save)				# W/m^2
            self.LW_flux_net			= np.zeros(self.nlev_save)				# W/m^2
            self.LW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev_save])	# W/m^2/(band)
            self.LW_spectral_flux_down 	= np.zeros([self.nbands,self.nlev_save])	# W/m^2/(band)
            self.LW_heating				= np.zeros(self.nlev_save)				# K/day
            self.SW_flux_up 			= np.zeros(self.nlev_save)				# W/m^2
            self.SW_flux_down 			= np.zeros(self.nlev_save)				# W/m^2
            self.SW_flux_net			= np.zeros(self.nlev_save)				# W/m^2
            self.SW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev_save])	# W/m^2/(band)
            self.SW_spectral_flux_down 	= np.zeros([self.nbands,self.nlev_save])	# W/m^2/(band)
            self.SW_heating				= np.zeros(self.nlev_save)				# K/day
            self.flux_up_total			= np.zeros(self.nlev_save)				# W/m^2
            self.flux_down_total		= np.zeros(self.nlev_save)				# W/m^2
            self.net_flux				= np.zeros(self.nlev_save)				# W/m^2
            self.net_spectral_flux	 	= np.zeros([self.nbands,self.nlev_save])	# W/m^2/(band)
            self.net_heating 			= np.zeros(self.nlev_save) 				# K/day from socrates
            self.heat                   = np.zeros(self.nlev_save)               # K/day from *
            self.cff					= np.zeros([self.nbands,self.nlev_save]) # W/m^2/(band). Compute by recompiling SOCRATES after setting the flags on lines 346 and 347 of src/aux/l_run_cdf.F90 to TRUE.
            self.LW_flux_up_i 			= np.zeros([self.nbands,self.nlev_save]) # W/m^2/(band)

        # Cloud flags (socrates/bin/rad_pcf.f90) and input arrays
        self.do_cloud = do_cloud
        if self.do_cloud :            
            self.cloud_scheme 			= 2	# -C 2 = ip_cloud_mix_max : maximum/random overlap in a mixed column 
            self.cloud_representation   = 1 # -K 1 = ip_cloud_homogen : ice and water mixed homogeneously 
            self.droplet_type           = 5 # -d 5 
            self.solver                 = 16 # -v 16 = i_solver : solver used for the two-stream calculations, chosen from those defined in solver_pcf.f90.
            # The below three variables should be zero until a cloud scheme forms clouds
            self.re                     = np.zeros(self.nlev_save) # Effective radius of the droplets [m]
            self.lwm                    = np.zeros(self.nlev_save) # Liquid water mass fraction [kg/kg]
            self.clfr                   = np.zeros(self.nlev_save) # Water cloud fraction
            # Floats defined in SocRadConv and used in the cloud scheme
            self.effective_radius       = re
            self.liquid_water_fraction  = lwm
            self.cloud_fraction         = clfr
        else:
            self.cloud_scheme 			= 5 # -C 5 = ip_cloud_off : clear sky. Other flags will be ignored. In principle we shouldn't need to define the variables below in this case.
            self.cloud_representation   = 1 
            self.droplet_type           = 5  
            self.solver                 = 13
            self.re                     = np.zeros(self.nlev_save)
            self.lwm                    = np.zeros(self.nlev_save)
            self.clfr                   = np.zeros(self.nlev_save) 
            self.effective_radius       = 0.0
            self.liquid_water_fraction  = 0.0
            self.cloud_fraction         = 0.0

        if self.alpha_cloud < 1.0e-20:
            self.effective_radius       = 0.0
            self.liquid_water_fraction  = 0.0
            self.cloud_fraction         = 0.0

    #New contructor based on toml file
    @classmethod
    def from_file(atm, file: str, band_edges: list, vol_mixing: dict = {}, vol_partial: dict = {}):

        with open(file, 'r') as f:
          cfg = toml.load(f)

        # Set parameters according to toml file if key present
        # Otherwise set default values
        T_surf     = cfg['atmos']['T_surf']     if "T_surf"     in cfg['atmos'] else 0.0
        P_surf     = cfg['atmos']['P_surf']     if "P_surf"     in cfg['atmos'] else 1.0e5
        P_top      = cfg['atmos']['P_top']      if "P_top"      in cfg['atmos'] else 1.0
        pl_radius  = cfg['planet']['pl_radius'] if "pl_radius"  in cfg['planet'] else 6.371e6
        pl_mass    = cfg['planet']['pl_mass']   if "pl_mass"    in cfg['planet'] else 5.972e24
        #optional arguments
        req_levels = cfg['atmos']['req_levels'] if "req_levels" in cfg['atmos'] else 100
        alpha      = cfg['atmos']['alpha']      if "alpha"      in cfg['atmos'] else 0.
        trppT      = cfg['atmos']['trppT']      if "trppT"      in cfg['atmos'] else 290.0
        do_cloud   = cfg['atmos']['do_cloud']   if "do_cloud"   in cfg['atmos'] else False
        re         = cfg['atmos']['re']         if "re"         in cfg['atmos'] else 0.
        lwm        = cfg['atmos']['lwm']        if "lwm"        in cfg['atmos'] else 0.
        clfr       = cfg['atmos']['clfr']       if "clfr"       in cfg['atmos'] else 0.
        albedo_pl  = cfg['atmos']['albedo_pl']  if "albedo_pl"  in cfg['atmos'] else 0.175
        albedo_s   = cfg['atmos']['albedo_s']   if "albedo_s"   in cfg['atmos'] else 0.
        zenith_angle = cfg['atmos']['zenith_angle'] if "zenith_angle" in cfg['atmos'] else 54.74

        return atm(T_surf, P_surf, P_top, pl_radius, pl_mass,
                   band_edges, vol_mixing, vol_partial,
                   #optional arguments
                   req_levels=req_levels,
                   alpha_cloud=alpha,
                   trppT=trppT,
                   do_cloud=do_cloud,
                   re=re,
                   lwm=lwm,
                   clfr=clfr,
                   albedo_pl=albedo_pl,
                   albedo_s=albedo_s,
                   zenith_angle=zenith_angle)

    def setSurfaceTemperature(self, Tsurf: float):

        self.ts = Tsurf
        self.tmp[0] = Tsurf
        self.tmpl[0] = Tsurf

    def setSurfacePressure(self, Psurf: float):

        self.ps = Psurf
        self.p[0] = Psurf
        self.pl[0] = Psurf

    def setPlanetProperties(self, pl_radius:float, pl_mass:float):

        self.planet_radius = pl_radius
        self.planet_mass = pl_mass
        self.grav_s = phys.G*self.planet_mass/(self.planet_radius**2) # m s-2
        self.grav_z[0] = self.grav_s

    def setVolatiles(self, vol_mixing: dict):

        tot_mixing =  float(sum(vol_mixing.values()))  # Ensure mixing ratios add up to unity
        self.vol_list = {}
        for vol in vol_mixing.keys():
          self.vol_list[vol] = vol_mixing[vol]/tot_mixing

        # H2O floor to prevent NaNs
        self.vol_list["H2O"] = np.max( [ self.vol_list["H2O"], 1e-20 ] )

        # Update volatile surface partial pressure
        for vol in self.vol_list.keys():
            self.p_vol[vol][0] = self.ps * self.vol_list[vol]

    def setTropopauseTemperature(self):

        T_eqm = (self.instellation * self.inst_sf * (1.0 - self.albedo_pl) /phys.sigma)**(1.0/4.0)
        self.trppT = T_eqm * (0.5**0.25)  # radiative skin temperature

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


    def write_ncdf(self, fpath:str):
        """Write atmosphere arrays to a netCDF file

        Parameters
        ----------
            fpath : string
                Output filename
        """ 

        # ----------------------
        # Calculate gravity and height (in case it hasn't been done already)
        self.z = AtmosphericHeight(self, self.planet_mass, self.planet_radius)

        self.zl = np.zeros(len(self.z)+1)
        for i in range(1,len(self.z)):
            self.zl[i] = 0.5 * (self.z[i-1] + self.z[i])
        self.zl[0] = 2*self.z[0] - self.zl[1] # estimate TOA height

        # ----------------------
        # Prepare NetCDF
        if os.path.exists(fpath):
            os.remove(fpath)

        ds = nc.Dataset(fpath, 'w', format='NETCDF4')
        ds.description        = 'JANUS atmosphere data'
        ds.hostname           = str(platform.node())
        try:
            # Try to get the login using os.getlogin()
            username = os.getlogin()
        except OSError:
            # If os.getlogin() fails, try an alternative method
            username = pwd.getpwuid(os.getuid()).pw_name
        ds.username = str(username)
        ds.JANUS_version     = "0.1"
        ds.SOCRATES_version   = "2306"
        ds.platform           = str(platform.system())

        nlev_c = len(self.p)
        nlev_l = nlev_c + 1

        gas_list = [str(gas) for gas in self.vol_list.keys()]
        ngases = len(gas_list)

        nchars = 16

        # ----------------------
        # Create dimensions
        ds.createDimension('nlev_c', nlev_c)    # Cell centres
        ds.createDimension('nlev_l', nlev_l)    # Cell edges
        ds.createDimension('ngases', ngases)    # Gases
        ds.createDimension('nchars', nchars)    # Length of string containing gas names
        ds.createDimension('nbands', self.nbands) # Number of bands in the spectral file
        ds.createDimension('nedges', self.nbands+1) 

        # ----------------------
        # Scalar quantities  
        #    Create variables
        var_tstar   = ds.createVariable('tstar',         'f4');  var_tstar.units = "K"     # BOA LW BC
        var_ps      = ds.createVariable('ps',            'f4');  var_ps.units = "Pa"
        var_ptop    = ds.createVariable('ptop',          'f4');  var_ptop.units = "Pa" 
        var_trppP   = ds.createVariable('trppP',         'f4');  var_trppP.units = "Pa" 
        var_inst    = ds.createVariable("instellation",  'f4');  var_inst.units = "W m-2"  # Solar flux at TOA
        var_s0fact  = ds.createVariable("inst_factor",   'f4');                            # Scale factor applied to instellation
        var_znth    = ds.createVariable("zenith_angle",  'f4');  var_znth.units = "deg"    # Scale factor applied to instellation
        var_albbond = ds.createVariable("bond_albedo",   'f4');                            # Bond albedo used to scale-down instellation
        var_toah    = ds.createVariable("toa_heating",   'f4');  var_toah.units = "W m-2"  # TOA SW BC
        var_tmagma  = ds.createVariable("tmagma",        'f4');  var_tmagma.units = "K"    # Magma temperature
        var_tmin    = ds.createVariable("tfloor",        'f4');  var_tmin.units = "K"      # Minimum temperature
        var_tmax    = ds.createVariable("tceiling",      'f4');  var_tmax.units = "K"      # Maximum temperature
        var_plrad   = ds.createVariable("planet_radius", 'f4');  var_plrad.units = "m"     # Value taken for planet radius
        var_gsurf   = ds.createVariable("surf_gravity",  'f4');  var_gsurf.units = "m s-2" # Surface gravity
        var_albsurf = ds.createVariable("surf_albedo",   'f4');                            # Surface albedo
        var_sknd    = ds.createVariable("cond_skin_d"   ,'f4');  var_sknd.units = "m"      # Conductive skin thickness
        var_sknk    = ds.createVariable("cond_skin_k"   ,'f4');  var_sknk.units = "W m-1 K-1"    # Conductive skin thermal conductivity

        #     Store data
        var_tstar.assignValue(self.ts)
        var_inst.assignValue(self.instellation)  
        var_toah.assignValue(self.toa_heating)
        var_ps.assignValue(self.ps)
        var_ptop.assignValue(self.ptop)
        var_trppP.assignValue(self.trppP)
        var_znth.assignValue(self.zenith_angle)
        var_s0fact.assignValue(self.inst_sf)
        var_albbond.assignValue(self.albedo_pl)
        var_tmagma.assignValue(self.tmp_magma)
        var_tmin.assignValue(self.minT)
        var_tmax.assignValue(self.maxT)
        var_plrad.assignValue(self.planet_radius)
        var_gsurf.assignValue(self.grav_s)
        var_albsurf.assignValue(self.albedo_s)
        var_sknd.assignValue(self.skin_d)
        var_sknk.assignValue(self.skin_k)


        # ----------------------
        # Layer quantities  
        #    Create variables
        var_p     = ds.createVariable('p',       'f4', dimensions=('nlev_c'));           var_p.units = "Pa"
        var_pl    = ds.createVariable('pl',      'f4', dimensions=('nlev_l'));           var_pl.units = "Pa"
        var_tmp   = ds.createVariable('tmp',     'f4', dimensions=('nlev_c'));           var_tmp.units = "K"
        var_tmpl  = ds.createVariable('tmpl',    'f4', dimensions=('nlev_l'));           var_tmpl.units = "K"
        var_z     = ds.createVariable('z',       'f4', dimensions=('nlev_c'));           var_z.units = "m"
        var_zl    = ds.createVariable('zl',      'f4', dimensions=('nlev_l'));           var_zl.units = "m"
        var_grav  = ds.createVariable('gravity', 'f4', dimensions=('nlev_c'));           var_grav.units = "m s-2"
        var_cp    = ds.createVariable('cp',      'f4', dimensions=('nlev_c'));           var_cp.units = "J/(kg K)"
        var_xd    = ds.createVariable('xd',      'f4', dimensions=('nlev_c'));           var_xd.units = "none"
        var_xv    = ds.createVariable('xv',      'f4', dimensions=('nlev_c'));           var_xv.units = "none"

        var_mmw   = ds.createVariable('mmw',     'f4', dimensions=('nlev_c'));           var_mmw.units = "kg mol-1"
        var_gases = ds.createVariable('gases',   'S1', dimensions=('ngases', 'nchars'))  # Names of gases
        var_mr    = ds.createVariable('x_gas',   'f4', dimensions=('nlev_c', 'ngases'))  # Mixing ratios per level
        var_cmr   = ds.createVariable('x_cond',  'f4', dimensions=('nlev_c', 'ngases'))  # Condensate mixing ratios per level
        var_pvol  = ds.createVariable('p_vol',   'f4', dimensions=('nlev_c', 'ngases')); var_pvol.units = "Pa"  # Gas phase partial pressures
        var_plvol = ds.createVariable('pl_vol',  'f4', dimensions=('nlev_l', 'ngases')); var_pvol.units = "Pa"
        if self.bands_set:
            var_fdl   = ds.createVariable('fl_D_LW', 'f4', dimensions=('nlev_l'));           var_fdl.units = "W m-2"
            var_ful   = ds.createVariable('fl_U_LW', 'f4', dimensions=('nlev_l'));           var_ful.units = "W m-2"
            var_fnl   = ds.createVariable('fl_N_LW', 'f4', dimensions=('nlev_l'));           var_fnl.units = "W m-2"
            var_fds   = ds.createVariable('fl_D_SW', 'f4', dimensions=('nlev_l'));           var_fds.units = "W m-2"
            var_fus   = ds.createVariable('fl_U_SW', 'f4', dimensions=('nlev_l'));           var_fus.units = "W m-2"
            var_fns   = ds.createVariable('fl_N_SW', 'f4', dimensions=('nlev_l'));           var_fns.units = "W m-2"
            var_fd    = ds.createVariable('fl_D',    'f4', dimensions=('nlev_l'));           var_fd.units = "W m-2"
            var_fu    = ds.createVariable('fl_U',    'f4', dimensions=('nlev_l'));           var_fu.units = "W m-2"
            var_fn    = ds.createVariable('fl_N',    'f4', dimensions=('nlev_l'));           var_fn.units = "W m-2"
            var_hr    = ds.createVariable('rad_hr',  'f4', dimensions=('nlev_c'));           var_hr.units = "K day-1"
            var_sful  = ds.createVariable('Sfl_U_LW','f4', dimensions=('nbands', 'nlev_l')); var_sful.units = "W m-2"
            var_sfus  = ds.createVariable('Sfl_U_SW','f4', dimensions=('nbands', 'nlev_l')); var_sfus.units = "W m-2"
            var_sfdl  = ds.createVariable('Sfl_D_LW','f4', dimensions=('nbands', 'nlev_l')); var_sfdl.units = "W m-2"
            var_sfds  = ds.createVariable('Sfl_D_SW','f4', dimensions=('nbands', 'nlev_l')); var_sfds.units = "W m-2"
            var_edge  = ds.createVariable('band_edges','f4', dimensions=('nedges'));         var_edge.units = "nm"

        if self.has_contfunc:
            var_cff   = ds.createVariable('cff',     'f4', dimensions=('nbands', 'nlev_c')); var_cff.units = "W m-2 m-1" ## units??

        var_re    = ds.createVariable('re',      'f4', dimensions=('nlev_c'));           var_re.units = "m"
        var_lwm   = ds.createVariable('lwm',     'f4', dimensions=('nlev_c'));           var_lwm.units = "kg kg-1"
        var_clfr  = ds.createVariable('clfr',    'f4', dimensions=('nlev_c'));           var_clfr.units = "none"

        #     Store data
        var_p[:]    = self.p[:]
        var_pl[:]   = self.pl[:]
        var_tmp[:]  = self.tmp[:]
        var_tmpl[:] = self.tmpl[:]
        var_z[:]    = self.z[:]
        var_zl[:]   = self.zl[:]
        var_mmw[:]  = self.mu[:]
        var_grav[:] = self.grav_z[:]
        var_cp[:]   = self.cp[:]
        var_xd[:]   = self.xd[:]
        var_xv[:]   = self.xv[:]

        var_gases[:] = np.array([ [c for c in gas.ljust(nchars)[:nchars]] for gas in gas_list ] , dtype='S1')
        var_mr[:]    = np.array([ [ self.x_gas[gas][i] for i in range(nlev_c-1,-1,-1) ] for gas in gas_list  ]).T
        var_cmr[:]   = np.array([ [ self.x_cond[gas][i] for i in range(nlev_c-1,-1,-1) ] for gas in gas_list  ]).T
        var_pvol[:]  = np.array([ [ self.p_vol[gas][i] for i in range(nlev_c-1,-1,-1) ] for gas in gas_list  ]).T
        var_plvol[:] = np.array([ [ self.pl_vol[gas][i] for i in range(nlev_l-1,-1,-1) ] for gas in gas_list  ]).T

        if  self.bands_set:
            var_fdl[:] = self.LW_flux_down[:]
            var_ful[:] = self.LW_flux_up[:]
            var_fnl[:] = self.LW_flux_net[:]

            var_fds[:] = self.SW_flux_down[:]
            var_fus[:] = self.SW_flux_up[:]
            var_fns[:] = self.SW_flux_net[:]

            var_fd[:] = self.flux_down_total[:]
            var_fu[:] = self.flux_up_total[:]
            var_fn[:] = self.net_flux[:]

            var_hr[:] = self.net_heating[:]

            var_sful[:,:]  = self.LW_spectral_flux_up[:,:]
            var_sfus[:,:]  = self.SW_spectral_flux_up[:,:]
            var_sfdl[:,:]  = self.LW_spectral_flux_down[:,:]
            var_sfds[:,:]  = self.SW_spectral_flux_down[:,:]

            var_edge[:]    = self.band_edges[:]

        if self.has_contfunc:
            var_cff[:,:]   = self.cff[:,:]

        var_re[:]   = self.re[:]
        var_lwm[:]  = self.lwm[:]
        var_clfr[:] = self.clfr[:]

        # ----------------------
        # Close
        ds.close()


