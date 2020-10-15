import numpy as np
import thermo
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

import standard_fluid_config as std
import injectors


class Pintle():
	def __init__(self, oxidiser_injector, fuel_injector, n_oxidiser_holes):
		"""Tool for sizing pintle injector uisng a pintle tip with discrete orifices and oxidiser centric injection

		:param oxidiser_injector: Oxidiser injector object 
		:param fuel_injector: Fuel injector object
		:param n_oxidiser_holes: number of orifices in the pintle tip
		"""		
		self.oxidiser_injector = oxidiser_injector
		self.fuel_injector = fuel_injector
		self.n_oxidiser_holes = n_oxidiser_holes
		self.oxidiser_injector.injector()                  # compute injector properties 
		self.fuel_injector.injector()                      # compute injector properties 

	def momentum_ratio(self):
		axial_momentum = self.fuel_injector.massflow * self.fuel_injector.velocity 
		
		radial_momentum = self.oxidiser_injector.massflow * self.oxidiser_injector.velocity * self.n_oxidiser_holes
	
		self.tmr = radial_momentum/axial_momentum
		self.efficiency = (-2*self.tmr**2 + 1.4*self.tmr + 92)/100

	def iterator(self, tmr_range):
		'''
		optimises pressure drops to get within given TMR range
		Increses pressure from the standard pressure drop 
		'''

		self.momentum_ratio() 								# initiate first momentum ratio calc
 
		while self.tmr > max(tmr_range) or self.tmr < min(tmr_range): 
			if self.tmr > max(tmr_range):
				self.fuel_injector.pressuredrop += 0.005e5
				self.fuel_injector.injector()

			if self.tmr < min(tmr_range):
				self.oxidiser_injector.pressuredrop += 0.005e5
				self.oxidiser_injector.injector()

			self.momentum_ratio()

class IjectorThermal():
	def __init__(self, thermal_conductivity, max_wall_thickness, fluid_temperature, fluid_pressure, fluid_massflow, fluid, fluid_mixture, hydrolic_diameter, massflow, chamber_pressure, velocity, gas_temperature=0):
		"""
		Estimates wall temperature at injection element. Using convective and radiative heat transfer 
		"""		
		self.thermal_conductivity = thermal_conductivity
		self.max_wall_thickness = max_wall_thickness
		self.fluid = thermo.Mixture(fluid, ws=fluid_mixture, T=fluid_temperature, P=fluid_pressure)
		self.flow_velocity = velocity
		self.chamber_diameter = std.chamber_diameter
		self.throat_diameter = std.throat_diameter
		self.chamber_pressure = chamber_pressure
		self.massflow = massflow
		self.hydrolic_diameter = hydrolic_diameter
		self.fluid_massflow = fluid_massflow
		self.gas_temperature = gas_temperature

		self.cea = std.cea

	def heat_trans_coeff_gas(self, mach, wall_temperature):
		gamma = self.cea.gamma
		t = self.cea.T_static
		if self.gas_temperature != 0:
			t= self.gas_temperature
		p = self.chamber_pressure
		mu = self.cea.mu
		cp = self.cea.Cp
		Pr = self.cea.Pr
		local_area = np.pi*self.chamber_diameter**2/4
		throat_area = np.pi*self.throat_diameter**2/4

		cstar = p * np.pi * self.throat_diameter**2 / 4 / self.massflow
		s = (0.5*wall_temperature/t * (1+(gamma-1)/2 * mach*mach) + 0.5)**(-0.68) * (1+(gamma-1)/2 * mach*mach)**(-0.12)
		halpha = 0.026/((self.throat_diameter)**0.2) * (mu**0.2)*cp/(Pr**0.6) * (p/cstar)**0.8 * (throat_area/local_area)**0.9 * s

		return halpha

	def adiabatic_wall_temp(self, mach):
		# Assumes turbulent flow in the chamber
		r = self.cea.Pr**(1/3) 					
		t_aw = self.cea.T_static * (1 + r*(self.cea.gamma - 1)/2 * mach**2) / (1 + (self.cea.gamma - 1)/2 * mach**2)

		if self.gas_temperature != 0:
			t_aw = self.gas_temperature
		
		return t_aw

	def radiation(self, mach):
		T_local = self.cea.T_static/(1 + (self.cea.gamma-1)/2 * mach**2)
		
		if self.gas_temperature != 0:
			T_local = self.gas_temperature

		p_co2 = self.cea.mole_fractions[1]['*CO2'][0] * self.chamber_pressure
		p_h2o = self.cea.mole_fractions[1]['H2O'][0] * self.chamber_pressure
		q_r_co2 = 4 * (p_co2/1e5*self.chamber_diameter/2)**0.3 * (T_local/100)**3.5
		q_r_h2o = 5.74 * (p_h2o/1e5*self.chamber_diameter/2)**0.3 * (T_local/100)**3.5

		return q_r_co2 + q_r_h2o

	def heat_trans_coeff_coolant(self, wall_temperature):
		flowvelocity = self.fluid_massflow/(self.fluid.rho * np.pi*(self.hydrolic_diameter/2)**2)
		Pr = self.fluid.Pr
		Re = self.fluid.rho*flowvelocity*hydrolic_diameter/self.fluid.mu
		k = self.fluid.Cp*self.fluid.mu/Pr
		
		Nu = 0.023*Re**0.8*Pr**0.4
		halpha = Nu*k/self.hydrolic_diameter

		return halpha, Re, Nu

	def wall_temperature(self, max_iter=1000, tol=1e-6):
		wall_temperature = 300
		iteration = 0
		difference_wall = 1

		R = 8314.5/self.cea.MW
		speed_of_sound = np.sqrt(self.cea.gamma*R*self.cea.T_static)
		mach = self.flow_velocity/speed_of_sound
		adiabatic_wall_temperature = self.adiabatic_wall_temp(mach)

		while difference_wall > tol:
			halpha = self.heat_trans_coeff_gas(mach, wall_temperature)
			halpha_c, Re, Nu = self.heat_trans_coeff_coolant(wall_temperature)
			radiation = self.radiation(mach)

			heat_flux = (adiabatic_wall_temperature - self.fluid.T + radiation/halpha) / (1/halpha + self.max_wall_thickness/self.thermal_conductivity + 1/halpha_c)
			new_wall_temp = - ((heat_flux - radiation)/halpha - adiabatic_wall_temperature)
			coolant_wall_temp = -heat_flux*self.max_wall_thickness/self.thermal_conductivity + new_wall_temp
			difference_wall = abs(new_wall_temp - wall_temperature)

			iteration += 1
			if iteration > max_iter:
				raise ValueError('Non-convergence, iteration number exceeded ', max_iter)

			wall_temperature = new_wall_temp

		self.max_wall_temperature = new_wall_temp
		self.heat_flux = heat_flux
		self.coolant_wall_temp = coolant_wall_temp

#Injector Parameters 
pressuredrop = std.pre_injection_pressure - std.chamber_pressure 		# [Pa]
fuel_inlet_angle = np.pi/2
n_holes = 30

l_hole = 2.75e-3							# [m]
annulus_length = 2e-3 						# [m]
d_pintle = 30e-3							# [m]

liq_inj = injectors.LiquidInjector(['O2'], [1], std.ox_temperature, std.pre_injection_pressure, l_hole, std.ox_massflow/n_holes, pressuredrop, fuel_inlet_angle)
an_inj = injectors.AnnulusInjector(std.fuel_composition, std.fuel_mass_fraction, std.fuel_injection_temperature, std.pre_injection_pressure, annulus_length,d_pintle, std.ox_massflow, pressuredrop)


# Pintle optimisation
tmr_range = [0.97,1]
pintle = Pintle(liq_inj, an_inj, n_holes)
pintle.iterator(tmr_range)
print('pintle injector TMR:', pintle.tmr)
print('oxidiser pressure drop: ', pintle.oxidiser_injector.pressuredrop/1e6, '[MPa]')
print('fuel pressure drop: ',pintle.fuel_injector.pressuredrop/1e6, '[MPa]')
print('oxidiser hole diameter:', pintle.oxidiser_injector.diameter*1000, '[mm]')
print('annulus width:', pintle.fuel_injector.diameter*1000, '[mm]')


# Injeector Thermals
thermal_conductivity_copper = 343			# [W/m/K]
thermal_conductivity_stainless = 27			# [W/m/K]
wall_thickness = 3e-3						# [m]
velocity = 100								# [m/s]
hydrolic_diameter = 12.25e-3				# [m]
manifold_hydrolic_diameter = 15e-3			# [m]
gas_temperature = 2674						# [K]

pintle_thermal = IjectorThermal(thermal_conductivity_copper, wall_thickness, std.ox_temperature, std.pre_injection_pressure, std.ox_massflow, std.ox_composition, [1], hydrolic_diameter, std.total_massflow, std.chamber_pressure, velocity)
pintle_thermal.wall_temperature()
print('maximum pintle tip temperature: ', pintle_thermal.max_wall_temperature, 'K')

faceplate_thermal = IjectorThermal(thermal_conductivity_stainless, wall_thickness, std.fuel_injection_temperature, std.pre_injection_pressure, std.fuel_massflow, std.fuel_composition, std.fuel_mass_fraction,manifold_hydrolic_diameter, std.total_massflow, std.chamber_pressure, velocity, gas_temperature)
faceplate_thermal.wall_temperature()
print('maximum face plate temperature: ', faceplate_thermal.max_wall_temperature, 'K')
print('coolant side face plate temperature: ', faceplate_thermal.coolant_wall_temp, 'K')