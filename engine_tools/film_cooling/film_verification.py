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
		y = geometry[:,1]
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
		y = geometry[:,1]
		t_aw = np.ndarray(len(y))

		for i in range(len(y)):
			if mach[i] >= 1:
				r = Pr**(1/2) 
			else:
				r = Pr**(1/3) 	
			t_aw[i] = self.t_s * (1 + r*(self.gamma - 1)/2 * mach[i]**2) / (1 + (self.gamma - 1)/2 * mach[i]**2)
		
		return t_aw

if __name__ == "__main__":
	
	geometry = np.genfromtxt('sparrow_50bar.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
	isnetropic = Isentropic(std.chamber_pressure, std.cea.T_static, std.cea.gamma)
	mach = isnetropic.mach(geometry)
	t_aw = isnetropic.adiabatic_wall_temp(mach, geometry, std.cea.Pr)