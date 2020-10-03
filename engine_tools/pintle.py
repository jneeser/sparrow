import numpy as np
import thermo
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

import standard_fluid_config as std
import injectors


class Pintle():
	def __init__(self, oxidiser_injector, fuel_injector, n_oxidiser_holes):
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



pressuredrop = std.pre_injection_pressure - std.chamber_pressure
fuel_inlet_angle = np.pi/2
n_holes = 30
l_hole = (15 - 12.25)*10**(-3)

ox_inlet_angle = np.pi/6
annulus_length =3e-3
d_pintle = 30e-3

liq_inj = injectors.LiquidInjector(['O2'], [1], std.ox_temperature, std.pre_injection_pressure, l_hole, std.ox_massflow/n_holes, pressuredrop, fuel_inlet_angle)
liq_inj.injector()
print(liq_inj.mu)
print(liq_inj.velocity)
print(liq_inj.diameter*1000)

an_inj = injectors.AnnulusInjector(std.fuel_composition, std.fuel_mass_fraction, std.fuel_temperature, std.pre_injection_pressure, annulus_length,d_pintle, std.ox_massflow, pressuredrop, ox_inlet_angle)
an_inj.injector()
print(an_inj.mu)
print(an_inj.velocity)
print(an_inj.diameter*1000)


tmr_range = [0.97,1]
pintle = Pintle(liq_inj, an_inj, n_holes)
pintle.iterator(tmr_range)
print(pintle.tmr)
print(pintle.efficiency)
print(pintle.oxidiser_injector.pressuredrop/1e6)
print(pintle.fuel_injector.pressuredrop/1e6)