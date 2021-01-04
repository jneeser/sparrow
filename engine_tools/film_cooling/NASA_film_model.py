import numpy as np
import scipy as sp
import thermo
from matplotlib import pyplot as plt
import rocketcea

from heat_transfer import CEA
import film_length as fl
import injectors as inj


class Isentropic():
	def __init__(self, total_pressure, total_temperature, gamma):
		self.p = total_pressure
		self.t = total_temperature
		self.gamma = gamma

	def mach(self, geometry):
		y = geometry[:,1][::-1]
		throat_diameter = 2*min(geometry[:,1])
		throat_area = np.pi*throat_diameter**2/4
		mach = np.ndarray(len(y))
		initial = 10

		for i in range(len(y)):
			if abs(y[i] - throat_diameter/2) < 1e-6:
				initial = 0.000001
			local_area = y[i]**2 * np.pi
			mach_number = lambda M:  1/(M*M) * (2/(self.gamma+1) * (1 + (self.gamma-1)/2*M*M))**((self.gamma+1)/(self.gamma-1)) - (local_area/throat_area)**2
			mach[i] = sp.optimize.fsolve(mach_number, initial)

		return mach

	def pressure(self,mach):
		return self.p/((1 + (self.gamma-1)/2 * mach**2)**(self.gamma/(self.gamma-1)))

	def temperature(self,mach):
		return self.t/(1 + (self.gamma-1)/2 * mach**2)

	def adiabatic_wall_temp(self, mach, geometry, Pr):
		# Assumes turbulent flow in the chamber and laminar flow after the throat
		T_aw = np.ndarray(len(geometry[:,1]))
		for i in range(len(T_aw)):
			if mach[i] >= 1:
				r = Pr**(1/2) 
			else:
				r = Pr**(1/3) 	
			T_aw[i] = self.t * (1 + r*(self.gamma - 1)/2 * mach[i]**2) / (1 + (self.gamma - 1)/2 * mach[i]**2)
		
		return T_aw


class FilmCooling():
	def __init__(self, coolant, cea, massflow, film_massflow, chamber_pressure, radius_at_injection, geometry):
		self.coolant = coolant 
		self.cea = cea
		self.massflow = massflow
		self.film_massflow = film_massflow
		self.chamber_pressure = chamber_pressure
		self.ri = radius_at_injection
		self.geometry = geometry

	def local_conditions(self, mach):
		self.u_e = mach*np.sqrt(self.cea.gamma*8314.472/self.cea.MW*self.cea.T_static)
		self.T_local = self.cea.T_static/(1 + (self.cea.gamma-1)/2 * mach**2)
		self.P_local = self.chamber_pressure/((1 + (self.cea.gamma-1)/2 * mach**2)**(self.cea.gamma/(self.cea.gamma-1)))
		self.rho_local = self.P_local / (8314.472/self.cea.MW * self.cea.T_static)
		self.T_if = self.coolant.Tc
		self.B = self.cea.Cp*(self.cea.T_static - self.T_if) / (self.coolant.Hvap_Tb - self.coolant.H)
	
	def liquid_lenght(self):
		
		if self.coolant.phase == 'g':
			raise ValueError('Fluid is not liquid')
		
		delta = 1 
		ri = self.ri*100/2.54								# conversion to [in]
		u_e = self.u_e / 0.3048								# conversion to [ft/s]
		film_massflow = self.film_massflow / 0.4536			# comversion to [lbm/s]
		rho = self.rho_local * 0.06242796					# conversion to [lbm/ft^3]
		
		sigma = self.coolant.SurfaceTension(T=self.coolant.T) * 0.0685217810		# conversion to [lbf/ft]
		Xe = delta*(rho/32.174)**0.5 * u_e * (self.cea.T_static/self.T_if)**0.25 / sigma
		Xr = Xe * sigma
		
		#from Liquid Rocket Engine Self-Cooled Combustion Chambers NASA 1977
		Xe_list = np.array([0,2,4,6,8,10,12,14,16,18,20])*1e4
		A_list = [0,0.1,0.2,0.3,0.38,0.45,0.5,0.56,0.6,0.66,0.7]
		par_list = [1,1.6,2.2,3.1,4,4.5,5,5.5,6,6.5,7]
		A_interp = sp.interpolate.interp1d(Xe_list,A_list,fill_value='extrapolate')
		par_interp = sp.interpolate.interp1d(Xe_list,par_list,fill_value='extrapolate')
		A = A_interp(Xe)			
		par = par_interp(Xe)
		a = par*(1+3*Xr**(-0.8))

		Re_gas = self.rho_local*self.u_e*2*self.ri / self.cea.mu
		St = 0.318*Re_gas**(-0.5)*self.cea.Pr**(-0.6)
		V = np.pi*rho*u_e/144 * 2*ri * St * self.B * a  		

		L = 1/A * np.log(1 + film_massflow * A / V) / 100 * 2.54 	# conversion to [m]

		return L 		

	def nasa_liquid_film(self, x_bar, phi_r=0.025):

		film_massflow = self.film_massflow / 0.4536			# comversion to [lbm/s]
		massflow = self.massflow / 0.4536					# comversion to [lbm/s]
		
		We_L = film_massflow*(1 / (0.6 * self.B/(1+self.B)) - 1)
		We_Wc = (massflow - film_massflow)/film_massflow * (2*phi_r*x_bar/self.ri) + We_L/film_massflow
			
		if massflow < (We_L+0.6*film_massflow):
			theta = 0.6 + 0.263 * (massflow-We_L)/film_massflow
		else:
			theta = 0.758

		eta = (theta*(1 + We_Wc*np.sqrt(1- We_L/(massflow-film_massflow)) - (phi_r*x_bar/self.ri)**2))**(-1)

		Cpv = self.coolant.Cpl
		h_total = self.cea.Cp*self.cea.T_static
		h_e = self.cea.Cp*self.T_local
		h_aw = h_total - eta*(h_total - self.coolant.H) - (1-self.cea.Pr**(1/3))*(h_total - h_e)
		h_c_sv = self.coolant.Hvap_Tb
		T_aw = (h_aw - eta*h_c_sv + eta*Cpv*self.T_if + (1-eta)*(self.cea.Cp*self.cea.T_static - h_e)) / (eta*Cpv + (1-eta)*self.cea.Cp)

		return T_aw

	def nasa_gaseous_film(self, x_bar, v_injector, d_injector):

		if self.coolant.phase == 'l':
			raise ValueError('Fluid is not gaseous')

		def f(uc_ue):
			#from Liquid Rocket Engine Self-Cooled Combustion Chambers NASA 1977
			if uc_ue < 1:
				f = (uc_ue)**1.5
			else: 
				f_list = [1,1.07,1.1,1.1,1.07,1.04,0.97,0.9,0.85,0.8,0.75,0.6,0.55,0.5,0.45]
				uc_ue_list = [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5,10]
				f_interp = sp.interpolate.interp1d(uc_ue_list,f_list,fill_value='extrapolate')
				f = f_interp(uc_ue)
			return f

		def eta_film(We_Wc):
			#from Liquid Rocket Engine Self-Cooled Combustion Chambers NASA 1977
			if We_Wc > 1.4:
				eta = 1.32 / (1+We_Wc)
			elif We_Wc > 0.06: 
				We_Wc_list = [0.06,0.2,0.4,0.8,1.2,1.4]
				eta_list = [1,0.9,0.8,0.7,0.6,0.55]
				eta_interp = sp.interpolate.interp1d(We_Wc_list,eta_list,fill_value='extrapolate')
				eta = eta_interp(We_Wc)
			else:
				eta = 1

			return eta

		uc_ue = v_injector/self.u_e
		psi_r = 0.1*uc_ue / ((self.coolant.rho/self.rho_local)**0.15 * (self.coolant.rho*v_injector*d_injector/self.coolant.mu)**0.25 * f(uc_ue))
		We_Wc = (self.massflow - self.film_massflow)/self.film_massflow * (2*psi_r*x_bar/(self.ri - d_injector) - (psi_r*x_bar/(self.ri - d_injector))**2)
		
		eta = eta_film(We_Wc)
		#print(eta)
		delta_h = self.cea.Cp * (self.cea.T_static - self.T_local)
		T_aw = self.cea.T_static - (eta*self.coolant.Cpg*(self.cea.T_static-self.coolant.T) + (1-self.cea.Pr**(1/3))*delta_h) / (eta*self.coolant.Cpg + (1-eta)*self.cea.Cp)

		return T_aw

	def T_aw(self, film_start, film_end, mach, T_aw_uncooled, n_holes, chamber_pressure):
		"""[summary]
		:param film_start: start index in 'geometry' of the film cooling
		:param film_end: expected index in 'geometry' of the film cooling for local geomoetry refinement 
		:param geometry: chamber geometry
		"""	
		geom_refinement = 10
		throat_idx = np.where(self.geometry[:,1]==min(self.geometry[:,1]))[0][0]*geom_refinement
		start = int(film_start*geom_refinement)
		end = int(film_end*geom_refinement)
		rg = fl.refine(self.geometry,geom_refinement)
		rx = rg[:,1]
		dx = rg[1,0]-rg[0,0]
		x = fl.x_bar(rx,start,dx,throat_idx)
		x_bar_arr = x.compute(start,end,start)

		pressure_drop = self.coolant.P - chamber_pressure
		liq_inj = inj.LiquidInjector(['c2h5oh', 'h2o'], [0.8,0.2], self.coolant.T, self.coolant.P, 0.6e-3, self.film_massflow/n_holes, pressure_drop, np.pi/6)
		liq_inj.injector()
		self.local_conditions(mach[0])
		L = self.liquid_lenght()

		points = np.arange(film_start, film_end, 1)
		T_aw_arr = np.ndarray(len(points))
		T_aw_cooled = T_aw_uncooled.copy()
		for p in points:
			i = p - film_start	
			self.local_conditions(mach[p])
			if self.geometry[i,0] > L:
				T_aw_arr[i] = T_aw_cooled[p]
			else: 
				T_aw_arr[i] = self.nasa_liquid_film(x_bar_arr[i*geom_refinement])

		T_aw_cooled[film_start:film_end] = T_aw_arr

		return T_aw_cooled



if __name__ == "__main__":

	OF = 1.61	
	oxidiser = 'LOX'
	ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])
	geometry = np.genfromtxt('sparrow_50bar.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
	massflow = 5.8
	film_massflow = 0.5
	chamber_pressure = 50e5
	mach = 0.2
	radius_at_injection = 60e-3

	cea = CEA(ethanol90, oxidiser, chamber_pressure)
	cea.metric_cea_output('throat', OF, 12)
	isnetropic = Isentropic(chamber_pressure, cea.T_static, cea.gamma)
	mach = isnetropic.mach(geometry)[::-1]
	T_aw_uncooled = isnetropic.adiabatic_wall_temp(mach, geometry, cea.Pr)

	coolant = thermo.Chemical('C2H5OH', P=60e5, T=350)

	film = FilmCooling(coolant, cea, massflow, film_massflow, chamber_pressure, geometry[42,1], geometry)
	#film.local_conditions(mach)
	#print(film.liquid_lenght())
	#print(film.nasa_liquid_film(0.3))
	#print(film.nasa_gaseous_film(0.2,10,0.4e-3))
	T_aw_cooled = film.T_aw(42, 70, mach,T_aw_uncooled, 40, chamber_pressure)

	f, axes = plt.subplots(4, 1)
	axes[0].plot(geometry[:,0], geometry[:,1]*1000)
	axes[0].set_ylabel('contour height [mm]')

	axes[1].plot(geometry[:,0], T_aw_uncooled)
	axes[1].set_ylabel('T_aw ')

	axes[2].plot(geometry[:,0], T_aw_cooled)
	axes[2].set_ylabel('T_aw cooled')

	axes[3].plot(geometry[:,0], mach)
	axes[3].set_ylabel('mach')

	plt.xlabel('x coordiante [m]')
	plt.show()