'''
##########################
For DARE internal use only 
##########################

AUTHOR: Jonthan Neeser
DATE: 26.08.2020
'''

import numpy as np
import thermo
from scipy.optimize import fsolve
from rocketcea.cea_obj import CEA_Obj

#INJECTOR CLASSES
#import injectors

class Pressurefed():
	#TODO add inulated piston and effect of pressurant temperature 
	#TODO add chemical gas generator (solid propellant)
	#TODO add compressibility effects on gas
	def __init__(self, oxidiser, fuel, pressurant, massflow_oxidiser, massflow_fuel, tankpressure, storagepressure, temperature_oxidiser, temperature_fuel, temperature_pressurant, burntime, piston=False):
		self.oxidiser = thermo.Chemical(oxidiser, T=temperature_oxidiser, P=tankpressure)
		self.fuel = thermo.Chemical(fuel,T=temperature_fuel, P=tankpressure)
		self.pressurant = thermo.Chemical(pressurant, T=temperature_pressurant, P=storagepressure)
		self.tankpressure = tankpressure
		self.pressurant_pressure = storagepressure
		self.massflow_fuel = massflow_fuel
		self.massflow_oxidiser = massflow_oxidiser
		self.burntime = burntime


	def pressurant_mass(self, margin = 0.05):
		volume = self.massflow_oxidiser*self.burntime/self.oxidiser.rho + self.massflow_fuel*self.burntime/self.fuel.rho
		volume *= (1+margin)
		gamma = self.pressurant.Cpg/self.pressurant.Cvg

		final_pressurant_pressure = self.oxidiser.P
		pressurant_mass = self.oxidiser.P*volume / (self.pressurant.R_specific*self.pressurant.T/2 * (1 + (final_pressurant_pressure/self.pressurant.P)**((gamma-1)/gamma))) / (1 - (final_pressurant_pressure/self.pressurant.P)**(1/gamma))
		pressurant_tank_volume = pressurant_mass*self.pressurant.R_specific*self.pressurant.T/self.pressurant.P

		avg_prop_temperature = (self.massflow_fuel*self.fuel.T + self.massflow_oxidiser*self.oxidiser.T)/(self.massflow_oxidiser+self.massflow_fuel)
		pressurant_mass_cooled = self.oxidiser.P*volume / (self.pressurant.R_specific*avg_prop_temperature/2 * (1 + (final_pressurant_pressure/self.pressurant.P)**((gamma-1)/gamma))) / (1 - (final_pressurant_pressure/self.pressurant.P)**(1/gamma))

		K = pressurant_mass_cooled/pressurant_mass

		return pressurant_mass_cooled, pressurant_tank_volume, K 

	def pressurant_mass_sutton(self, margin=0.05):
		volume = self.massflow_oxidiser*self.burntime/self.oxidiser.rho + self.massflow_fuel*self.burntime/self.fuel.rho
		volume *= (1+margin)
		gamma = self.pressurant.Cpg/self.pressurant.Cvg
		
		pressurant_mass = self.oxidiser.P * volume / (self.pressurant.R_specific * self.pressurant.T) * gamma / (1 - self.oxidiser.P/self.pressurant.P)

		avg_prop_temperature = (self.massflow_fuel*self.fuel.T + self.massflow_oxidiser*self.oxidiser.T)/(self.massflow_oxidiser+self.massflow_fuel)
		pressurant_mass_cooled = self.oxidiser.P * volume / (self.pressurant.R_specific * avg_prop_temperature) * gamma / (1 - self.oxidiser.P/self.pressurant.P)

		K = pressurant_mass_cooled/pressurant_mass

		pressurant_tank_volume = pressurant_mass*self.pressurant.R_specific*self.pressurant.T/self.pressurant.P

		return pressurant_mass_cooled, pressurant_tank_volume, K

	def tank_mass(self, rocket_diameter, ullage=0.05):
		# from NASA subsystem mass database 
		_, pressurant_tank_volume, _ = self.pressurant_mass_sutton()
		prop_volume = self.massflow_oxidiser*self.burntime/self.oxidiser.rho + self.massflow_fuel*self.burntime/self.fuel.rho
		prop_volume *= (1+ullage)

		m_prop_tanks = prop_volume*self.tankpressure/(6.43e4)
		m_pressurant_tank = pressurant_tank_volume*self.pressurant_pressure/(6.43e4)

		return m_prop_tanks, m_pressurant_tank



class Isentropic():
	def __init__(self, static_pressure, static_temperature, gamma):
		self.p_s = static_pressure
		self.t_s = static_temperature
		self.gamma = gamma

	def mach(self, chamber_area, throat_area, diverging):
		if diverging == 1:
			initial_guess = 10
		else:
			initial_guess = 0.0001
		mach = lambda M:  1/(M*M) * (2/(self.gamma+1) * (1 + (self.gamma-1)/2*M*M))**((self.gamma+1)/(self.gamma-1)) - (chamber_area/throat_area)**2
		mach_number = fsolve(mach, initial_guess)
		return mach_number

	def pressure(self,mach):
		return self.p_s/((1 + (self.gamma-1)/2 * mach**2)**(self.gamma/(self.gamma-1)))

	def temperature(self,mach):
		return self.t_s/(1 + (self.gamma-1)/2 * mach**2)

		

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
	def __init__(self, coolant, coolant_massfraction, coolant_massflow, total_massflow, fuel, oxidiser, mixture_ratio, chamber_pressure, coolant_temperature, coolant_pressure, geometry, wall_thickness, thermal_conductivity):
		"""[summary]
		Coolant flow properties stored as total conditions
		Hot gas properties from CEA, currently assumes constant gas properties in chamber 
		USE SI UNITS
		"""
		self.chamber_pressure = chamber_pressure
		self.mixture_ratio = mixture_ratio
		self.massflow = total_massflow
		self.coolant_massflow = coolant_massflow
		self.coolant = thermo.Mixture(coolant, ws = coolant_massfraction, T=coolant_temperature, P=coolant_pressure)
		self.chamber_diameter = 2*geometry[0,1]
		self.throat_diameter = 2*min(geometry[:,1])
		self.expansion_ratio = np.pi*geometry[-1][1]**2/(np.pi*self.throat_diameter**2/4)
		self.wall_thickness = wall_thickness
		self.thermal_conductivity = thermal_conductivity

		# get hot gas properties from CEA
		self.cea = CEA(fuel, oxidiser, self.chamber_pressure)
		self.cea.metric_cea_output('chamber', self.mixture_ratio, self.expansion_ratio)

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

	def adiabatic_wall_temp(self, mach, diverging):
		# Assumes turbulent flow in the chamber and laminar flow after the throat
		if diverging == 1:
			r = self.cea.Pr**(1/2) 
		else:
			r = self.cea.Pr**(1/3) 					
		t_aw = self.cea.T_static * (1 + r*(self.cea.gamma - 1)/2 * mach**2) / (1 + (self.cea.gamma - 1)/2 * mach**2)
		
		return t_aw

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

	def heat_trans_coeff_coolant(self, hydrolic_diameter, wall_temperature, x_coordinate, y_coordinate):
		coolant_area = hydrolic_diameter / 2 * y_coordinate * 2 * np.pi
		flowvelocity = self.coolant_massflow/(self.coolant.rho * coolant_area)
		Pr = self.coolant.Pr
		Re = self.coolant.rho*flowvelocity*hydrolic_diameter/self.coolant.mu
		k = self.coolant.Cp*self.coolant.mu/Pr

		Nu = 0.023*Re**0.8*Pr**0.4#*(self.coolant.T/wall_temperature) ** (0.57 - 1.59*hydrolic_diameter/x_coordinate)
		halpha = Nu*k/hydrolic_diameter

		#G = self.coolant_massflow/coolant_area
		#halpha = 0.029*self.coolant.Cp*self.coolant.mu**0.2/Pr**(2/3) * (G**0.8/hydrolic_diameter**0.2) * (self.coolant.T/wall_temperature) ** 0.55

		return halpha, Re, Nu

	def iterator(self, y_coordinate, x_coordinate, hydrolic_diameter, section_lenght, initial_guess, mach, adiabatic_wall_temperature ,max_iter=1000, tol=1e-6):
		wall_temperature = 300
		iteration = 0
		difference = 1

		while difference > tol:
			halpha = self.heat_trans_coeff_gas(mach, wall_temperature, y_coordinate)
			halpha_c, Re, Nu = self.heat_trans_coeff_coolant(hydrolic_diameter, wall_temperature, x_coordinate, y_coordinate)
			radiation = self.radiation(y_coordinate, mach)

			heat_flux = (adiabatic_wall_temperature - self.coolant.T + radiation/halpha) / (1/halpha + self.wall_thickness/self.thermal_conductivity + 1/halpha_c)
			new_wall_temp = - ((heat_flux - radiation)/halpha - adiabatic_wall_temperature)

			difference = abs(new_wall_temp - wall_temperature)
			iteration += 1
			if iteration > max_iter:
				raise ValueError('Non-convergence, iteration number exceeded ', max_iter)

			wall_temperature = new_wall_temp

		T_new = self.coolant.T + heat_flux*2*np.pi*y_coordinate*section_lenght / (self.coolant_massflow*self.coolant.Cp) 
		dp = self.pressure_drop(0, hydrolic_diameter, section_lenght, y_coordinate)
		self.coolant.calculate(P=self.coolant.P-dp, T=T_new)

		return heat_flux, wall_temperature, Re, Nu, radiation

	def heatflux(self, hydrolic_diameter, geometry, max_temperature=0, optimise=False):
		"""determines heat flux along the entire geometry starting from the nozzle end. Calls iterator function for all grid points. Only use for engine with radial cooling jacket. Can optimise cooling flow hydrolic diameter for a maximum wall temperature 

		:param hydrolic_diameter: hyrolic diamter of cooling passage
		:type hydrolic_diamter: array
		:param geometry: array of x and y coordinates of the chamber grid points
		:type geometry: array
		:param max_temperature: ONLY FOR OPTIMISATION, max allowable wall temp 
		:type max_temperature: float
		:param optimise: set to TRUE to optimise the non-cylindical cooling channel sections for a certian wall temperature 
		:type optimise: boolean

		###################################
		OUTPUTS (at each chamber location):
		###################################
		mach: 							local mach number 
		adiabatic_wall_temperature: 	adiabatic wall temperature 
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
		self.mach = np.ndarray(len(y))
		self.adiabatic_wall_temperature = np.ndarray(len(y))
		self.wall_temp = np.ndarray(len(y))
		self.coolant_temp = np.ndarray(len(y))
		self.coolant_pressure = np.ndarray(len(y))
		self.q = np.ndarray(len(y))
		self.q_rad = np.ndarray(len(y))
		self.coolant_Re = np.ndarray(len(y))
		self.coolant_Nu = np.ndarray(len(y))
		self.optimised_hydrolic_diameter = hydrolic_diameter
		self.T_chamber = np.ndarray(len(y))
		self.P_chamber = np.ndarray(len(y))


		# initial guess for mach-area relation
		initial_guess = np.ndarray(len(y))					
		diverging = True
		for i in range(len(y)):
			initial_guess[i] = diverging
			if abs(y[i] - self.throat_diameter/2) < 1e-6:
				diverging = False

		local = Isentropic(self.chamber_pressure, self.cea.T_static, self.cea.gamma)

		# Iterate over each chamber lcoation 
		for i in range(len(y)):
			if i == 0:
				section_length = 0
			else:
				section_length = np.sqrt((x[i] - x[i-1])**2 + (y[i]-y[i-1])**2)

			local_area = np.pi*y[i]**2
			mach = local.mach(local_area, np.pi*self.throat_diameter**2/4, initial_guess[i])
			adiabatic_wall_temperature = self.adiabatic_wall_temp(mach, initial_guess[i])
			
			self.mach[i] = mach
			self.adiabatic_wall_temperature[i] = adiabatic_wall_temperature
			self.P_chamber[i] = local.pressure(mach)
			self.T_chamber[i] = local.temperature(mach)
		
			q, wall_temp, Re, Nu, radiation = self.iterator(y[i], x[i], hydrolic_diameter[i], section_length, initial_guess[i], mach, adiabatic_wall_temperature)

			# if optimise = True optimise cooling jacket geometry
			if optimise and y[i] < self.chamber_diameter:
				while wall_temp > max_temperature:
					hydrolic_diameter[i] -= 0.05e-3
					q, wall_temp, Re, Nu, radiation = self.iterator(y[i], x[i], hydrolic_diameter[i], section_length, initial_guess[i], mach, adiabatic_wall_temperature)
			
			self.q[i] = q  
			self.q_rad[i] = radiation
			self.wall_temp[i] = wall_temp
			self.coolant_temp[i] = self.coolant.T    
			self.coolant_pressure[i] = self.coolant.P    
			self.coolant_Re[i] = Re
			self.coolant_Nu[i] = Nu
			self.optimised_hydrolic_diameter = hydrolic_diameter
		


	
