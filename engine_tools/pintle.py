import numpy as np
import thermo
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
		radial_momentum = self.oxidiser_injector.massflow * self.oxidiser_injector.velocity * self.n_oxidiser_holes
		axial_momentum = self.fuel_injector.massflow * self.fuel_injector.velocity 
		self.tmr = radial_momentum/axial_momentum

		radial_local_area = np.pi * self.oxidiser_injector.diameter**2/4
		axial_local_area = self.fuel_injector.diameter*self.oxidiser_injector.diameter	
		self.lmr = self.oxidiser_injector.fluid.rho*self.oxidiser_injector.velocity**2*radial_local_area / (self.fuel_injector.fluid.rho*self.fuel_injector.velocity**2*axial_local_area)
 		
		alpha = 0.7
		beta = 2.0
		self.spray_angle = alpha * np.arctan(beta * self.lmr)

	def pintle_optimiser(self, tmr_range):
		"""
		optimises pressure drops to get within given TMR range
		Increses pressure from the standard pressure drop value
		"""
		self.momentum_ratio() 								# initiate first momentum ratio calc
 
		while self.tmr > max(tmr_range) or self.tmr < min(tmr_range): 
			if self.tmr > max(tmr_range):
				self.fuel_injector.pressuredrop += 0.005e5
				self.fuel_injector.injector()

			if self.tmr < min(tmr_range):
				self.oxidiser_injector.pressuredrop += 0.005e5
				self.oxidiser_injector.injector()

			self.momentum_ratio()

	def pintle_sensitivity(self, uncertainty=0.1):
		ox_pressure_range = np.linspace(self.oxidiser_injector.pressuredrop*(1-uncertainty), self.oxidiser_injector.pressuredrop*(1+uncertainty), 100)
		fuel_pressure_range = np.linspace(self.fuel_injector.pressuredrop*(1-uncertainty), self.fuel_injector.pressuredrop*(1+uncertainty), 100)

		spray_angles = np.ndarray((len(ox_pressure_range), len(fuel_pressure_range)))
		tmr = np.ndarray((len(ox_pressure_range), len(fuel_pressure_range)))
		
		cd_ox = self.oxidiser_injector.mu
		cd_fuel = self.fuel_injector.mu
		
		for i in range(len(ox_pressure_range)):
			for j in range(len(fuel_pressure_range)):
				v_ox = cd_ox*np.sqrt(2*ox_pressure_range[i]/self.oxidiser_injector.fluid.rho)
				v_fuel = cd_fuel*np.sqrt(2*fuel_pressure_range[j]/self.fuel_injector.fluid.rho)
				m_ox = self.oxidiser_injector.fluid.rho * v_ox * self.n_oxidiser_holes * np.pi*self.oxidiser_injector.diameter**2/4
				m_fuel = self.fuel_injector.fluid.rho * v_fuel * np.pi*self.fuel_injector.D*self.fuel_injector.diameter

				radial_local_area = np.pi * self.oxidiser_injector.diameter**2/4
				axial_local_area = self.fuel_injector.diameter*self.oxidiser_injector.diameter	
				lmr = self.oxidiser_injector.fluid.rho*v_ox**2*radial_local_area / (self.fuel_injector.fluid.rho*v_fuel**2*axial_local_area)

				alpha = 0.7
				beta = 2.0
				spray_angles[i][j] = alpha * np.arctan(beta * lmr)

				radial_momentum = m_ox * v_ox 
				axial_momentum = m_fuel * v_fuel
				tmr[i][j] = radial_momentum/axial_momentum

		x,y = np.meshgrid(ox_pressure_range,fuel_pressure_range)

		fig, ax = plt.subplots()
		ax.set_title('spray angles [-]')
		c = ax.pcolormesh(x/1e6, y/1e6, np.degrees(spray_angles))
		plt.scatter(self.oxidiser_injector.pressuredrop/1e6, self.fuel_injector.pressuredrop/1e6, self.spray_angle, color='red', label='design point')
		ax.axis([np.min(x)/1e6, np.max(x)/1e6, np.min(y)/1e6, np.max(y)/1e6])
		plt.xlabel('dp oxidiser [MPa]')
		plt.ylabel('dp fuel [MPa]')
		plt.legend(loc='best')
		fig.colorbar(c)

		plt.show()

		fig, ax = plt.subplots()
		ax.set_title('TMR [-]')
		c = ax.pcolormesh(x/1e6, y/1e6, tmr)
		plt.scatter(self.oxidiser_injector.pressuredrop/1e6, self.fuel_injector.pressuredrop/1e6, self.tmr, color='red', label='design point')
		ax.axis([np.min(x)/1e6, np.max(x)/1e6, np.min(y)/1e6, np.max(y)/1e6])
		plt.xlabel('dp oxidiser [MPa]')
		plt.ylabel('dp fuel [MPa]')
		plt.legend(loc='best')
		fig.colorbar(c)

		plt.show()


#Injector Parameters 
pressuredrop = std.pre_injection_pressure - std.chamber_pressure 		# [Pa]
inlet_angle = np.pi/2

n_holes = 36
l_hole = 2.75e-3							# [m]
annulus_length = 2e-3 						# [m]
d_pintle = 30e-3							# [m]

liq_inj = injectors.LiquidInjector(['o2'], [1], std.ox_temperature, std.pre_injection_pressure, l_hole, std.ox_massflow/n_holes, pressuredrop, inlet_angle)
an_inj = injectors.AnnulusInjector(['c2h5oh', 'h2o'], [0.9,0.1], std.fuel_injection_temperature, std.pre_injection_pressure, annulus_length, d_pintle, std.fuel_massflow, pressuredrop)

# Pintle optimisation 
tmr_range = [0.9,1.1]
pintle = Pintle(liq_inj, an_inj, n_holes)
pintle.pintle_optimiser(tmr_range)
print('pintle injector TMR:', pintle.tmr)
print('pintle injector spray cone half angle:', np.degrees(pintle.spray_angle))
print('oxidiser pressure drop: ', pintle.oxidiser_injector.pressuredrop/1e6, '[MPa]')
print('fuel pressure drop: ',pintle.fuel_injector.pressuredrop/1e6, '[MPa]')
print('oxidiser hole diameter:', pintle.oxidiser_injector.diameter*1000, '[mm]')
print('annulus width:', pintle.fuel_injector.diameter*1000, '[mm]')

# pintle spray angle sensitivity 
pintle.pintle_sensitivity()


