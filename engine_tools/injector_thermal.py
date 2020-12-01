import numpy as np
import thermo
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

import standard_fluid_config as std

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

	def heat_trans_coeff_coolant(self, wall_temperature, flowvelocity):
		Pr = self.fluid.Pr
		Re = self.fluid.rho*flowvelocity*self.hydrolic_diameter/self.fluid.mu
		k = self.fluid.Cp*self.fluid.mu/Pr
		
		Nu = 0.023*Re**0.8*Pr**0.4
		halpha = Nu*k/self.hydrolic_diameter

		return halpha, Re, Nu

	def wall_temperature(self, flowvelocity, max_iter=1000, tol=1e-6):
		wall_temperature = 300
		iteration = 0
		difference_wall = 1

		R = 8314.5/self.cea.MW
		speed_of_sound = np.sqrt(self.cea.gamma*R*self.cea.T_static)
		mach = self.flow_velocity/speed_of_sound
		adiabatic_wall_temperature = self.adiabatic_wall_temp(mach)

		while difference_wall > tol:
			halpha = self.heat_trans_coeff_gas(mach, wall_temperature)
			halpha_c, Re, Nu = self.heat_trans_coeff_coolant(wall_temperature, flowvelocity)
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
		self.halpha_c = halpha_c
		self.halpha = halpha


if __name__ == "__main__":

	# Pintle Tip Thermals
	thermal_conductivity_copper = 343			# [W/m/K]
	wall_thickness = 3e-3						# [m]
	velocity = 100								# [m/s]
	hydrolic_diameter = 24.5e-3/2				# [m]
	flowvelocity = std.ox_massflow/(std.liquid_ox.rho * np.pi*(hydrolic_diameter/2)**2)

	pintle_thermal = IjectorThermal(thermal_conductivity_copper, wall_thickness, std.ox_temperature, std.pre_injection_pressure, std.ox_massflow, std.ox_composition, [1], hydrolic_diameter, std.total_massflow, std.chamber_pressure, velocity)
	pintle_thermal.wall_temperature(flowvelocity)
	print('maximum pintle tip temperature: ', pintle_thermal.max_wall_temperature, 'K')


	#Faceplate Thermals 
	thermal_conductivity_stainless = 28			# [W/m/K]
	wall_thickness = 3e-3						# [m]
	velocity = 100								# [m/s]
	gas_temperature = 2376						# [K]
	local_OF = 0.97

	chamber_radi = np.linspace(16e-3, 65e-3, 100)
	flow_height = 3e-3							# [m]
	width = 2*np.pi*chamber_radi
	areas = width*flow_height
	flow_velocity = std.fuel_massflow / std.liquid_fuel.rho / areas
	hydrolic_diameter = 2*flow_height*areas / (areas + flow_height)

	coolant_side_faceplate_temp = np.ndarray(len(chamber_radi))
	gas_side_faceplate_temp = np.ndarray(len(chamber_radi))
	halpha_c = np.ndarray(len(chamber_radi))

	for i in range(len(chamber_radi)):
		faceplate_thermal = IjectorThermal(thermal_conductivity_stainless, wall_thickness, std.fuel_injection_temperature, std.pre_injection_pressure, std.fuel_massflow, std.fuel_composition, std.fuel_mass_fraction,hydrolic_diameter[i], std.total_massflow, std.chamber_pressure, velocity)
		faceplate_thermal.wall_temperature(flow_velocity[i])
		coolant_side_faceplate_temp[i] = faceplate_thermal.coolant_wall_temp
		gas_side_faceplate_temp[i] = faceplate_thermal.max_wall_temperature
		halpha_c[i] = faceplate_thermal.halpha_c


	print('maximum face plate temperature: ', max(gas_side_faceplate_temp), 'K')
	print('maximum coolant side face plate temperature: ', max(coolant_side_faceplate_temp), 'K')

	plt.plot(chamber_radi*1e3, gas_side_faceplate_temp, label="maximum faceplate temperature")
	plt.plot(chamber_radi*1e3, coolant_side_faceplate_temp, label="coolant side faceplate temperature")
	plt.xlabel("chamber radius [mm]")
	plt.ylabel("temperature [K]")
	plt.legend(loc="best")
	plt.grid()
	plt.show()

	plt.plot(chamber_radi*1e3, flow_velocity)
	plt.xlabel("chamber radius [mm]")
	plt.ylabel("flow velocity [m/s]")
	plt.grid()
	plt.show()


	# export geometric parameters in mm for catia import 
	with open('coolant_halpha.csv', 'w', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(["chamber radius [m]", "coolant heat transfer coefficient [W/m^2/K]"]
		)
		for idx in range(len(chamber_radi)):
			writer.writerow([chamber_radi[idx], halpha_c[idx]])

	