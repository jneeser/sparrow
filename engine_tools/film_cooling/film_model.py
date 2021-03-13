import numpy as np
import scipy as sp
import thermo
from matplotlib import pyplot as plt
import rocketcea

from heat_transfer import CEA
import injectors as inj


class FilmCooling():
	def __init__(self, coolant, cea, massflow, chamber_pressure, coolant_velocity):
		self.coolant = coolant
		self.cea = cea
		self.massflow = massflow
		self.chamber_pressure = chamber_pressure
		self.coolant_velocity = coolant_velocity

	def local_conditions(self, mach):
		self.ue = mach*np.sqrt(self.cea.gamma*8314.472/self.cea.MW*self.cea.T_static)
		self.T_local = self.cea.T_static/(1 + (self.cea.gamma-1)/2 * mach**2)
		self.P_local = self.chamber_pressure/((1 + (self.cea.gamma-1)/2 * mach**2)**(self.cea.gamma/(self.cea.gamma-1)))
		self.rho_local = self.P_local / (8314.472/self.cea.MW * self.cea.T_static)

	def liquid_film_cooling(self, chamber_radius, q_conv, q_rad, film_massflow):
		Tb = sum([self.coolant.Tbs[i]*self.coolant.ws[i] for i in range(len(self.coolant.ws))])
		hfg = thermo.phase_change.Clapeyron(self.coolant.T, self.coolant.Tc, self.coolant.Pc) * self.coolant.MW 		# conversion to J/kg
		hfg = hfg + (self.coolant.Tc - self.coolant.T) * self.coolant.Cpl
		mv = (q_conv+q_rad) / hfg

		Gmean = self.rho_local*self.ue * (self.ue - self.coolant_velocity)/self.ue
		Re = Gmean * 2 * chamber_radius / self.cea.mu

		func = lambda f: 1/np.sqrt(f) - 1.93*np.log10(Re*np.sqrt(f)) + 0.537
		f = sp.optimize.fsolve(func, 0.00001)
		St0 = f/2 / (1.2 + 11.8 * np.sqrt(f/2) * (self.cea.Pr - 1) * self.cea.Pr**(-1/3))
		Kt = 1.4
		h0 = Gmean * self.cea.Cp * St0 * Kt
		
		F = mv / Gmean
		Mg_Mc = self.cea.MW / self.coolant.MW

		F_St0 = self.cea.Cp/hfg * (self.T_local - self.coolant.Tc + q_rad/h0)
		St_St0 = np.log(1 + F_St0 * (Mg_Mc)**0.6) / (F_St0 * (Mg_Mc)**0.6)

		self.eta = St_St0
		self.film_length = film_massflow/(2*np.pi*chamber_radius)/mv
		self.mv = mv


	def film_effectiveness(self, film_start, section_length, mach, q_conv, q_rad, film_massflow, geometry):
		"""[summary]

		:param film_start: 		start index of film cooling on chamber geometry
		:param section_length	array of section lenghts along contour 
		:param mach: 			array of local mach numbers
		:param q_conv: 			array of convective heat flux, starting at injector
		:param q_rad: 			array of radiation heat flux, starting at injector
		:param film_massflow:	initial coolant mass flow 
		"""	

		self.eta_arr = np.ndarray(len(geometry[:,1]))
		y = geometry[:,1]

		for i in range(len(y)):
			self.local_conditions(mach[i])
			if film_massflow < 0:
				self.eta_arr[i] = 1
			elif i >= film_start:
				self.liquid_film_cooling(y[i], q_conv[i], q_rad[i], film_massflow)
				self.eta_arr[i] = self.eta
				film_massflow -= self.mv * 2*np.pi*y[i] * section_length[i]
				print(self.film_length*1000)
				print(q_conv[i]/1e6)
			else:
				self.eta_arr[i] = 1	


if __name__ == "__main__":
	

	OF = 1.61	
	oxidiser = 'LOX'
	ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])
	geometry = np.genfromtxt('sparrow_50bar.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
	massflow = 5.8
	film_massflow = 0.4878
	chamber_pressure = 50e5
	mach = 0.2
	coolant_velocity = 70
	q_rad = 2e6
	q_conv = 8.1e6


	cea = CEA(ethanol90, 'LOX', chamber_pressure)
	cea.metric_cea_output('throat', OF, 12)

	coolant = thermo.Mixture(['C2H5OH','H2O'], ws=[0.9,0.1], P=60e5, T=288)

	film = FilmCooling(coolant, cea, massflow, chamber_pressure, coolant_velocity)
	film.local_conditions(mach)
	film.liquid_film_cooling(0.0607, q_conv, q_rad, film_massflow)
	print(film.film_length)

	mach = np.ones(100)*0.2
	q_rad = np.ones(100)*1.8e6
	q_conv = np.ones(100)*15e6
	section_length = np.ones(100)*0.01
	y = np.ones(100)*0.065
	x = np.arange(0,1,0.01)

	#film.film_effectiveness(50, section_length, mach, q_conv, q_rad, film_massflow, y)


	#plt.plot(x, film.eta_arr)
	#plt.show()