'''
##########################
For DARE internal use only 
##########################

AUTHOR: Jonthan Neeser
DATE: 13.10.2020
'''

import numpy as np
import thermo
from scipy.optimize import fsolve
from rocketcea.cea_obj import CEA_Obj


class Isentropic():
	def __init__(self, static_pressure, static_temperature, gamma):
		self.p_s = static_pressure
		self.t_s = static_temperature
		self.gamma = gamma

	def mach(self, geometry):
		y = geometry[:,1][::-1]
		throat_diameter = 2*min(geometry[:,1])
		throat_area = np.pi*throat_diameter**2/4
		mach = np.ndarray(len(y))
		initial_guess = np.ndarray(len(y))	
		
		guess = 10
		for i in range(len(y)):
			initial_guess[i] = guess
			if abs(y[i] - throat_diameter/2) < 1e-6:
				guess = 0.000001

		for i in range(len(y)):
			local_area = y[i]**2 * np.pi
			mach_number = lambda M:  1/(M*M) * (2/(self.gamma+1) * (1 + (self.gamma-1)/2*M*M))**((self.gamma+1)/(self.gamma-1)) - (local_area/throat_area)**2
			mach[i] = fsolve(mach_number, initial_guess[i])

		return mach

	def pressure(self,mach):
		return self.p_s/((1 + (self.gamma-1)/2 * mach**2)**(self.gamma/(self.gamma-1)))

	def temperature(self,mach):
		return self.t_s/(1 + (self.gamma-1)/2 * mach**2)

	def adiabatic_wall_temp(self, mach, geometry, Pr):
		# Assumes turbulent flow in the chamber and laminar flow after the throat
		y = geometry[:,1][::-1]
		t_aw = np.ndarray(len(y))

		for i in range(len(y)):
			if mach[i] >= 1:
				r = Pr**(1/2) 
			else:
				r = Pr**(1/3) 	
			t_aw[i] = self.t_s * (1 + r*(self.gamma - 1)/2 * mach[i]**2) / (1 + (self.gamma - 1)/2 * mach[i]**2)
		
		return t_aw


class CEA():
	def __init__(self, fuel, oxidiser, chamber_pressure):
		"""[summary]
		Calcualtes hot gas properties using CEA and converts to SI units. If errors occur with FORTRAN, resart and try again!
		:param chamber_pressure in [Pa]
		:type fuel = string
		:type oxidiser = string
		"""
		self.chamber_pressure = chamber_pressure
		self.imperial_pressure = 0.000145038*self.chamber_pressure 			#conversion to psia
		self.ispObj = CEA_Obj( oxName=oxidiser, fuelName=fuel)

	def chamber_gas_properties(self, mixture_ratio, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Chamber_Transport(Pc=self.imperial_pressure, MR=mixture_ratio, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_Chamber_MolWt_gamma(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)

	def throat_gas_properties(self, mixture_ratio, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Throat_Transport(Pc=self.imperial_pressure, MR=mixture_ratio, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_Throat_MolWt_gamma(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)
	
	def exit_gas_properties(self, mixture_ratio, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Exit_Transport(Pc=self.imperial_pressure, MR=mixture_ratio, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_exit_MolWt_gamma(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)


	def metric_cea_output(self, location, mixture_ratio, expansion_ratio):
		if location == 'chamber':
			self.chamber_gas_properties(mixture_ratio, expansion_ratio)
		elif location == 'throat':
			self.throat_gas_properties(mixture_ratio, expansion_ratio)
		elif location == 'exit':
			self.exit_gas_properties(mixture_ratio, expansion_ratio)
		else:
			raise ValueError('Invalid location, use "chamber," "throat" or "exit"')

		self.isp, self.cstar, _ =  self.ispObj.getFrozen_IvacCstrTc(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)
		
		self.cstar = self.cstar * 0.3048																										# coversion to m/s
		self.mole_fractions = self.ispObj.get_SpeciesMoleFractions(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio, frozen=0, frozenAtThroat=0, min_fraction=5e-05)
		self.Cp = self.Cp * 4186.8																												# coversion to J/gk/K
		self.mu = self.visc * 0.0001																											# coversion to Pa*s
		self.k = self.cond * 418.4e-3																											# coversion to W/m/K
		self.T_static = self.ispObj.get_Tcomb(Pc=self.imperial_pressure, MR=mixture_ratio)*0.555556										        # coversion to K		
		


class Heattransfer():
	#TODO add curvature correction factors
	#TODO add support for only cooled chamber
	#TODO add varying gas properties in chamber
	#TODO add support for thermal barrier coating
	def __init__(self, coolant, coolant_massfraction, coolant_massflow, total_massflow, cea, chamber_pressure, coolant_temperature, coolant_pressure, geometry, number_of_channels, thermal_conductivity, k_tbc=0, t_tbc=0):
		"""[summary]
		Coolant flow properties stored as total conditions
		Hot gas properties from CEA, currently assumes constant gas properties in chamber 
		"""
		self.chamber_pressure = chamber_pressure
		self.massflow = total_massflow
		self.coolant_massflow = coolant_massflow
		self.coolant = thermo.Mixture(coolant, ws = coolant_massfraction, T=coolant_temperature, P=coolant_pressure)
		self.chamber_diameter = 2*geometry[0,1]
		self.throat_diameter = 2*min(geometry[:,1])
		self.expansion_ratio = np.pi*geometry[-1][1]**2/(np.pi*self.throat_diameter**2/4)
		self.thermal_conductivity = thermal_conductivity
		self.coolant_species = coolant
		self.coolant_massfraction = coolant_massfraction
		self.number_of_channels = number_of_channels
		self.k_tbc = k_tbc
		self.t_tbc = t_tbc

		# get hot gas properties from CEA
		self.cea = cea

	def heat_trans_coeff_gas(self, mach, wall_temperature, y_coordinate):
		gamma = self.cea.gamma
		t = self.cea.T_static
		p = self.chamber_pressure
		mu = self.cea.mu
		cp = self.cea.Cp
		Pr = self.cea.Pr
		local_area = np.pi*y_coordinate**2
		throat_area = np.pi*self.throat_diameter**2/4

		cstar = p * np.pi * self.throat_diameter**2 / 4 / self.massflow
		s = (0.5*wall_temperature/t * (1+(gamma-1)/2 * mach*mach) + 0.5)**(-0.68) * (1+(gamma-1)/2 * mach*mach)**(-0.12)
		halpha = 0.026/((self.throat_diameter)**0.2) * (mu**0.2)*cp/(Pr**0.6) * (p/cstar)**0.8 * (throat_area/local_area)**0.9 * s

		return halpha

	def pressure_drop(self, surface_roughness, hydrolic_diameter, section_lenght, y_coordinate):
		coolant_area = hydrolic_diameter / 2 * y_coordinate * 2 * np.pi
		flowvelocity = self.coolant_massflow/(self.coolant.rho * coolant_area)
		Re = self.coolant.rho*flowvelocity*hydrolic_diameter/self.coolant.mu
		fd = fsolve(lambda f: 1/(np.sqrt(f)) + 2*np.log10(surface_roughness/(3.7*hydrolic_diameter) + 2.51/(Re*np.sqrt(f))), 0.0000001)
		dp = fd*section_lenght/hydrolic_diameter*0.5*self.coolant.rho*flowvelocity**2 
		
		return dp
		
	def radiation(self, y_coordinate, mach):
		T_local = self.cea.T_static/(1 + (self.cea.gamma-1)/2 * mach**2)
		p_co2 = self.cea.mole_fractions[1]['*CO2'][0] * self.chamber_pressure
		p_h2o = self.cea.mole_fractions[1]['H2O'][0] * self.chamber_pressure
		q_r_co2 = 4 * (p_co2/1e5*y_coordinate)**0.3 * (T_local/100)**3.5
		q_r_h2o = 5.74 * (p_h2o/1e5*y_coordinate)**0.3 * (T_local/100)**3.5

		return q_r_co2 + q_r_h2o

	def heat_trans_coeff_coolant(self, hydrolic_diameter, wall_temperature):
		coolant_area = hydrolic_diameter**2/4 * np.pi * self.number_of_channels 
		flowvelocity = self.coolant_massflow/(self.coolant.rho * coolant_area)
		Pr = self.coolant.Pr
		Re = self.coolant.rho*flowvelocity*hydrolic_diameter/self.coolant.mu
		k = self.coolant.Cp*self.coolant.mu/Pr
		
		Nu = 0.023*Re**0.8*Pr**0.4
		halpha = Nu*k/hydrolic_diameter

		return halpha, Re, Nu

	def iterator(self, y_coordinate, hydrolic_diameter, section_length, wall_thickness, mach, adiabatic_wall_temperature ,max_iter=1000, tol=1e-6):
		wall_temperature = 600
		tbc_wall_temperature = 600
		iteration = 0
		difference_wall = 1

		while difference_wall > tol:
			halpha = self.heat_trans_coeff_gas(mach, tbc_wall_temperature, y_coordinate)
			halpha_c, Re, Nu = self.heat_trans_coeff_coolant(hydrolic_diameter, wall_temperature)
			radiation = self.radiation(y_coordinate, mach)

			if self.k_tbc == 0:
				# no thermal barrier coating 
				heat_flux = (adiabatic_wall_temperature - self.coolant.T + radiation/halpha) / (1/halpha + wall_thickness/self.thermal_conductivity + 1/halpha_c)
				new_wall_temp = - ((heat_flux - radiation)/halpha - adiabatic_wall_temperature)
				tbc_wall_temperature = new_wall_temp
			else:
				# with thermal barrier coating 
				heat_flux = (adiabatic_wall_temperature - self.coolant.T + radiation/halpha) / (1/halpha + wall_thickness/self.thermal_conductivity + 1/halpha_c + self.t_tbc/self.k_tbc)
				new_wall_temp = ((halpha*halpha_c*self.k_tbc*wall_thickness*adiabatic_wall_temperature + 
								self.coolant.T*halpha*halpha_c*self.thermal_conductivity*self.t_tbc + 
								halpha*self.thermal_conductivity*self.k_tbc*adiabatic_wall_temperature +
								self.coolant.T*halpha_c*self.thermal_conductivity*self.k_tbc) / 
								(halpha*halpha_c*self.thermal_conductivity*self.t_tbc + 
								halpha*halpha_c*self.k_tbc*wall_thickness + 
								halpha*self.thermal_conductivity*self.k_tbc + 
								halpha_c*self.thermal_conductivity*self.k_tbc)
				)
				tbc_wall_temp = - ((heat_flux - radiation)/halpha - adiabatic_wall_temperature)

			difference_wall = abs(new_wall_temp - wall_temperature)

			iteration += 1
			if iteration > max_iter:
				raise ValueError('Non-convergence, iteration number exceeded ', max_iter)

			wall_temperature = new_wall_temp

		return heat_flux, wall_temperature, tbc_wall_temp, Re, Nu, radiation, halpha

	def heatflux(self, hydrolic_diameter, geometry, wall_thickness, mach_numbers, adiabatic_wall_temperatures):
		"""determines heat flux along the entire geometry starting from the nozzle end. Calls iterator function for all grid points. Only use for engine with radial cooling jacket. Can optimise cooling flow hydrolic diameter for a maximum wall temperature 

		:param hydrolic_diameter: hyrolic diamter of cooling passage
		:type hydrolic_diamter: array
		:param geometry: array of x and y coordinates of the chamber grid points
		:type geometry: array
		:param wall_thickness: array of local inner wall thicknesses 
		:type wall_thickness: array
		:param mach_numbers: mach numbers along each geometry section starting from the nozzle
		:type mach_numbers: array
		:param adiabatic_wall_temperatures: adiabatic wall temperatures along each geometry section starting from the nozzle
		:type adiabatic_wall_temperatures: array


		###################################
		OUTPUTS (at each chamber location):
		###################################
		wall_temp:						wall temperature of the inner chamber 
		coolant_temp: 					bulk temperature of the coolant 
		coolant_pressure: 				pressure of the coolant 
		q: 								total heat exchanged betwen chamber and coolant 
		q_rad:							total radiative heat flux 
		coolant_Re:						Bulk Reynolds number in the cooling passage 
		coolant_Nu:						Bulk Nusselt number in the cooling passage 
		optimised_hydrolic_diameter:	Hydrolic dimaeter after optimisation 
		"""        
		y = geometry[:,1][::-1]
		x = geometry[:,0][::-1]

		# create empty output arrays
		self.wall_temp = np.ndarray(len(y))
		self.coolant_temp = np.ndarray(len(y))
		self.coolant_pressure = np.ndarray(len(y))
		self.q = np.ndarray(len(y))
		self.q_rad = np.ndarray(len(y))
		self.coolant_Re = np.ndarray(len(y))
		self.coolant_Nu = np.ndarray(len(y))
		self.halpha_gas = np.ndarray(len(y))
		self.tbc_wall_temp = np.ndarray(len(y))


		# Iterate over each chamber lcoation 
		for i in range(len(y)):
			if i == 0:
				section_length = 0
			else:
				section_length = np.sqrt((x[i] - x[i-1])**2 + (y[i]-y[i-1])**2)
		
			q, wall_temp, tbc_wall_temp, Re, Nu, radiation, halpha = self.iterator(y[i], hydrolic_diameter[i], section_length, wall_thickness[i], mach_numbers[i], adiabatic_wall_temperatures[i])
			T_new = self.coolant.T + q*2*np.pi*y[i]*section_length / (self.coolant_massflow*self.coolant.Cp) 
			self.coolant.calculate(P=self.coolant.P, T=T_new)

			self.q[i] = q  
			self.q_rad[i] = radiation
			self.wall_temp[i] = wall_temp
			self.coolant_temp[i] = self.coolant.T    
			self.coolant_pressure[i] = self.coolant.P    
			self.coolant_Re[i] = Re
			self.coolant_Nu[i] = Nu
			self.halpha_gas[i] = halpha
			self.tbc_wall_temp[i] = tbc_wall_temp
		


	
