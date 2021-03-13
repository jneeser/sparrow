'''
##########################
For DARE internal use only 
##########################

AUTHOR: Jonthan Neeser
DATE: 13.12.2020
'''

import numpy as np
import thermo
import scipy.optimize 
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
		mach_number = scipy.optimize.fsolve(mach, initial_guess)
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
	def __init__(self, coolant, coolant_massfraction, coolant_massflow, total_massflow, fuel, oxidiser, mixture_ratio, chamber_pressure, coolant_temperature, coolant_pressure, geometry, cooling_geometry, thermal_conductivity, method='standard-bartz', T_aw_cooled=0):
		"""[summary]
		Coolant flow properties stored as total conditions
		Hot gas properties from CEA, currently assumes constant gas properties in chamber 
		USE SI UNITS
		"""
		self.geometry = geometry
		self.cooling_geometry = cooling_geometry
		self.chamber_pressure = chamber_pressure
		self.mixture_ratio = mixture_ratio
		self.massflow = total_massflow
		self.coolant_massflow = coolant_massflow
		self.coolant = thermo.Mixture(coolant, ws = coolant_massfraction, T=coolant_temperature, P=coolant_pressure)
		self.chamber_diameter = 2*geometry[0,1]
		self.throat_diameter = 2*min(geometry[:,1])
		self.expansion_ratio = np.pi*geometry[-1][1]**2/(np.pi*self.throat_diameter**2/4)
		self.thermal_conductivity = thermal_conductivity
		self.coolant_species = coolant
		self.coolant_massfraction = coolant_massfraction
		self.method = method
		self.T_aw_cooled = T_aw_cooled

		# get hot gas properties from CEA
		self.cea = CEA(fuel, oxidiser, self.chamber_pressure)
		self.cea.metric_cea_output('throat', self.mixture_ratio, self.expansion_ratio)

	def heat_trans_coeff_gas(self, mach, wall_temperature, t_aw, y_coordinate):
		gamma = self.cea.gamma
		p = self.chamber_pressure
		mu = self.cea.mu
		cp = self.cea.Cp
		Pr = self.cea.Pr
		T_local = self.cea.T_static/(1 + (self.cea.gamma-1)/2 * mach**2)
		local_area = np.pi*y_coordinate**2
		throat_area = np.pi*self.throat_diameter**2/4
		
		cstar = p * np.pi * self.throat_diameter**2 / 4 / self.massflow

		if self.method == 'standard-bartz':
			s = (0.5*wall_temperature/self.cea.T_static * (1+(gamma-1)/2 * mach*mach) + 0.5)**(-0.68) * (1+(gamma-1)/2 * mach*mach)**(-0.12)
			halpha = 0.026/((self.throat_diameter)**0.2) * (mu**0.2)*cp/(Pr**0.6) * (p/cstar)**0.8 * (throat_area/local_area)**0.9 * s

		elif self.method == 'modified-bartz':
			T_f = 0.5*wall_temperature + 0.28*T_local + 0.22*t_aw
			G = self.massflow/local_area
			halpha = 0.026 * G**0.8/(self.throat_diameter)**0.2 * mu**0.2*cp/Pr**0.6 * (self.cea.T_static/T_f)**0.68

		elif self.method == 'cinjarew':
			eta_combustion = 0.92
			T_hg = T_local + 0.8*(self.cea.T_static*eta_combustion**2 - T_local)
			halpha = 0.01975 * self.cea.k**0.18*(self.massflow*cp)**0.82 / (2*y_coordinate)**1.82 * (T_hg/wall_temperature)**0.35

		else:
			raise ValueError('Invalid heat transfer method. Select: "standard-bartz", "modified-bartz" or "cinjarew"')

		return halpha

	def adiabatic_wall_temp(self, mach, diverging):
		# Assumes turbulent flow in the chamber and laminar flow after the throat
		if diverging == 1:
			r = self.cea.Pr**(1/2) 
		else:
			r = self.cea.Pr**(1/3) 					
		t_aw = self.cea.T_static * (1 + r*(self.cea.gamma - 1)/2 * mach**2) / (1 + (self.cea.gamma - 1)/2 * mach**2)
		
		return t_aw

	def pressure_drop(self, surface_roughness, section_length, section_number, y_coordinate):
		d_h = self.cooling_geometry.dhi_arr[section_number]
		A = self.cooling_geometry.Ai_arr[section_number]

		flowvelocity = self.coolant_massflow/(self.coolant.rho * A * 2 * self.cooling_geometry.N)
		Re = self.coolant.rho*flowvelocity*d_h/self.coolant.mu
		fd = scipy.optimize.fsolve(lambda f: -2*np.log10(surface_roughness/(3.7*d_h)+2.51/(Re*np.sqrt(f)))-1/np.sqrt(f), 0.0000001)
		dp = fd*section_length/d_h*0.5*self.coolant.rho*flowvelocity**2 
		
		return dp
		
	def radiation(self, y_coordinate, mach):
		T_local = self.cea.T_static/(1 + (self.cea.gamma-1)/2 * mach**2)
		P_local = self.chamber_pressure/((1 + (self.cea.gamma-1)/2 * mach**2)**(self.cea.gamma/(self.cea.gamma-1)))
		p_co2 = self.cea.mole_fractions[1]['*CO2'][0] * P_local
		p_h2o = self.cea.mole_fractions[1]['H2O'][0] * P_local
		q_r_co2 = 4 * (p_co2/1e5*y_coordinate)**0.3 * (T_local/100)**3.5
		q_r_h2o = 5.74 * (p_h2o/1e5*y_coordinate)**0.3 * (T_local/100)**3.5

		return q_r_h2o + q_r_co2

	def heat_trans_coeff_coolant(self, wall_temperature, coolant_wall_temperature, x_coordinate, y_coordinate, section_length, section_number):
		d_h = self.cooling_geometry.dhi_arr[section_number]
		A = self.cooling_geometry.Ai_arr[section_number]
		
		flowvelocity = self.coolant_massflow/(self.coolant.rho * A * 2 * self.cooling_geometry.N)
		Pr = 0.75 + 1.63/np.log(1+self.coolant.Pr/0.0015) 			# turbulent Pr correction

		Re = self.coolant.rho*flowvelocity*d_h/self.coolant.mu
		k = self.coolant.Cp*self.coolant.mu/Pr

		wall_fluid = thermo.Mixture(self.coolant_species, ws=self.coolant_massfraction, P=self.coolant.P, T=coolant_wall_temperature)
		Nu = 0.0208*Re**0.8*Pr**0.4*(1+0.01457*wall_fluid.mu/self.coolant.mu)  #Hess & Kunz relationship
		halpha = Nu * k / d_h
		
		#halpha = 0.023*self.coolant.Cp**0.333*k**0.667 / (self.coolant.mu**0.467*d_h**0.2) * (self.coolant_massflow/(np.pi/4 * d_h**2))**0.8  # McAdams
		#halpha = 0.023*k/d_h * (self.coolant.Cp / (k*self.coolant.mu))**0.4 * (self.coolant_massflow*d_h/A)**0.8
		#Nu = halpha / k * d_h

		return halpha, Re, Nu, flowvelocity

	def iterator(self, y_coordinate, x_coordinate, section_length, section_number, initial_guess, mach, t_aw ,eta=1, max_iter=1000, tol=1e-6):
		wall_temperature = 300
		coolant_wall_temperature = 300
		wall_thickness = self.cooling_geometry.wt1_arr[section_number]
		iteration = 0
		difference = 1

		while difference > tol:
			halpha = self.heat_trans_coeff_gas(mach, wall_temperature, t_aw, y_coordinate) * eta						# film cooling correction
			halpha_c, Re, Nu, flowvelocity = self.heat_trans_coeff_coolant(wall_temperature, coolant_wall_temperature, x_coordinate, y_coordinate, section_length, section_number)
			radiation = self.radiation(y_coordinate, mach)

			heat_flux = (t_aw - self.coolant.T + radiation/halpha) / (1/halpha + wall_thickness/self.thermal_conductivity + 1/halpha_c)
			new_wall_temp = - ((heat_flux - radiation)/halpha - t_aw)

			difference = abs(new_wall_temp - wall_temperature)

			iteration += 1
			if iteration > max_iter:
				raise ValueError('Non-convergence, iteration number exceeded ', max_iter)

			wall_temperature = new_wall_temp
			coolant_wall_temperature = -heat_flux*wall_thickness/self.thermal_conductivity + new_wall_temp

		T_new = self.coolant.T + heat_flux*2*np.pi*y_coordinate*section_length / (self.coolant_massflow*self.coolant.Cp) 
		dp = self.pressure_drop(6e-6, section_length, section_number, y_coordinate)
		self.coolant.calculate(P=self.coolant.P-dp, T=T_new)

		return heat_flux, wall_temperature, Re, Nu, flowvelocity, radiation, halpha, halpha_c

	def heatflux(self, geometry, eta_film=1, x_film_injection=0):
		"""
		OUTPUTS (at each chamber location):
		mach: 							local mach number 
		t_aw: 							adiabatic wall temperature 
		wall_temp:						wall temperature of the inner chamber 
		coolant_temp: 					bulk temperature of the coolant 
		coolant_pressure: 				pressure of the coolant 
		q: 								total heat exchanged betwen chamber and coolant 
		q_rad:							total radiative heat flux 
		coolant_Re:						Bulk Reynolds number in the cooling passage 
		coolant_Nu:						Bulk Nusselt number in the cooling passage 
		tbc_wall_temp:					Wall temoperature outside of thermal barrier coating 
		flowvelocity: 					Velocity of flow in the cooling channels 				
		"""        
		y = self.geometry[:,1][::-1]
		x = self.geometry[:,0][::-1]

		# create empty output arrays
		self.mach = np.ndarray(len(y))
		self.t_aw = np.ndarray(len(y))
		self.wall_temp = np.ndarray(len(y))
		self.coolant_temp = np.ndarray(len(y))
		self.coolant_pressure = np.ndarray(len(y))
		self.q = np.ndarray(len(y))
		self.q_rad = np.ndarray(len(y))
		self.coolant_Re = np.ndarray(len(y))
		self.coolant_Nu = np.ndarray(len(y))
		self.T_chamber = np.ndarray(len(y))
		self.P_chamber = np.ndarray(len(y))
		self.halpha_gas = np.ndarray(len(y))
		self.flowvelocity = np.ndarray(len(y))
		self.section_length = np.ndarray(len(y))
		self.halpha_coolant = np.ndarray(len(y))

		# initial guess for mach-area relation
		initial_guess = np.ndarray(len(y))					
		diverging = True
		for i in range(len(y)):
			initial_guess[i] = diverging
			if abs(y[i] - self.throat_diameter/2) < 1e-6:
				diverging = False

		local = Isentropic(self.chamber_pressure, self.cea.T_static, self.cea.gamma)

		for i in range(len(y)):
			if i == 0:
				section_length = 0
			else:
				section_length = np.sqrt((x[i] - x[i-1])**2 + (y[i]-y[i-1])**2)

			local_area = np.pi*y[i]**2
			mach = local.mach(local_area, np.pi*self.throat_diameter**2/4, initial_guess[i])
			t_aw = self.adiabatic_wall_temp(mach, initial_guess[i])
			eta = 1

			if type(self.T_aw_cooled) != int:
				t_aw = self.T_aw_cooled[i]

			if type(eta_film) != int:
				eta = eta_film[i]

			self.mach[i] = mach
			self.t_aw[i] = t_aw
			self.P_chamber[i] = local.pressure(mach)
			self.T_chamber[i] = local.temperature(mach)
			self.section_length[i] = section_length
		
			q, wall_temp, Re, Nu, flowvelocity, radiation, halpha, halpha_c = self.iterator(y[i], x[i], section_length, i, initial_guess[i], mach, t_aw, eta)

			self.q[i] = q  
			self.q_rad[i] = radiation
			self.wall_temp[i] = wall_temp
			self.coolant_temp[i] = self.coolant.T    
			self.coolant_pressure[i] = self.coolant.P    
			self.coolant_Re[i] = Re
			self.coolant_Nu[i] = Nu
			self.halpha_gas[i] = halpha
			self.flowvelocity[i] = flowvelocity
			self.halpha_coolant[i] = halpha_c



	
