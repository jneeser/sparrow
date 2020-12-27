import numpy as np
import thermo
import rocketcea
from matplotlib import pyplot as plt
import csv

import heat_transfer as ht
import geom_class as gc
import NASA_film_model as film

geometry = np.genfromtxt('trotti_contour.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
t = 2e-3
wt1 = 1e-3
wt2 = 0.01e-3
rf1 = 0.01e-3
rf2 = 0.01e-3
number_of_channels = 42
cooling_channels = gc.parameters(ri=25e-3,t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=number_of_channels)
cooling_channels.cooling_geometry(geometry[:,1][::-1])


thermal_conductivity = 392							# copper 
method = 'cinjarew'

# global properties
chamber_pressure = 40e5 			# [Pa]
fuel_inlet_pressure = 80e5			# [Pa]
fuel_temperature = 288				# [K]
ox_temperature = 288 				# [K]
expansion_ratio = 1.5

# CEA input values 
OF = 3.22	
oxidiser = 'GOX'
fuel = 'Kerosene'
cea = ht.CEA(fuel, oxidiser, chamber_pressure)
cea.metric_cea_output('chamber', OF, expansion_ratio)

# Thermo input values 
total_massflow = 0.508						# [kg/s]
fuel_massflow = 3
fuel_composition = ['H2O']
fuel_mass_fraction = [1]

# Film Cooling BETA
mu = 0.15
film_massflow = mu*total_massflow
isnetropic = film.Isentropic(chamber_pressure, cea.T_static, cea.gamma)
mach = isnetropic.mach(geometry)[::-1]
T_aw_uncooled = isnetropic.adiabatic_wall_temp(mach, geometry, cea.Pr)
coolant = thermo.Chemical('C2H5OH', P=chamber_pressure+10e5, T=288)
film = film.FilmCooling(coolant, cea, total_massflow, film_massflow, chamber_pressure, geometry[52,1], geometry)
T_aw_cooled = film.T_aw(film_start=52, film_end=80, mach=mach, T_aw_uncooled=T_aw_uncooled, n_holes=number_of_channels, chamber_pressure=chamber_pressure)[::-1]


heat = ht.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, fuel, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, cooling_channels, thermal_conductivity, method, T_aw_cooled)

heat.heatflux(geometry)

plt.rcParams.update({'font.size': 12})
f, axes = plt.subplots(4, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.t_aw)
axes[3].set_ylabel('adiabatic wall temperature [K]')

plt.xlabel('x coordiante [m]')
plt.show()

print(max(heat.wall_temp))

''' 
with open('heat_transfer_coefficients.csv', 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(["x coordinate","y_coordinate","gas heat transfer coefficient"])
	for i in range(len(geometry[:,1])):
		writer.writerow([geometry[i,0], geometry[i,1], heat.halpha_gas[len(geometry[:,1]) - i - 1]])
'''

plt.plot(geometry[:,0][::-1], heat.q/1e6, label='total heat flux')
#plt.plot(geometry[:,0][::-1], heat.halpha_gas*(heat.t_aw-heat.wall_temp)/1e6, label='convective heat flux')
#plt.plot(geometry[:,0][::-1], heat.q_rad/1e6, label='radiation')
plt.scatter([0.05, 0.15, 0.25, 0.385],np.array([1.289e7,1.545e7, 7.468e6,1.561e7])/1e6, color='red', label='experimental heat fluxes')
plt.grid()
plt.xlabel('x coordinate [m]')
plt.ylabel('heat flux [MW/m^2]')
plt.legend(loc='best')
plt.show()
