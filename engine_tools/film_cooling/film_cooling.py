import numpy as np
import thermo
import rocketcea
from matplotlib import pyplot as plt
import csv

import heat_transfer as ht
import geom_class as gc
import NASA_film_model as film
import film_model as fm

geometry = np.genfromtxt('sparrow_08atm.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
t = 2.4e-3
wt1 = 0.55e-3
wt2 = 0.6e-3
rf1 = 0.3e-3
rf2 = 0.4e-3
number_of_channels = 42
cooling_channels = gc.parameters(ri=25e-3,t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=number_of_channels)
cooling_channels.cooling_geometry(geometry[:,1][::-1])


thermal_conductivity = 21.2							# Inconel at 1000 K
method = 'cinjarew'

# global properties
chamber_pressure = 50e5 			# [Pa]
fuel_inlet_pressure = 70e5			# [Pa]
fuel_temperature = 288				# [K]
ox_temperature = 90 				# [K]
expansion_ratio = 10

# CEA input values 
OF = 1.49	
oxidiser = 'LOX'
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[80,20])  # new fuel blend for CEA
cea = ht.CEA(ethanol90, oxidiser, chamber_pressure)
cea.metric_cea_output('throat', OF, expansion_ratio)

# Thermo input values 
total_massflow = 5.8 						# [kg/s]
fuel_massflow = total_massflow / (1+OF)
mu = 0.15
film_massflow = fuel_massflow * mu / (1-mu)
print(film_massflow)

ox_massflow = total_massflow - fuel_massflow
fuel_massflow += film_massflow
fuel_composition = ['C2H5OH']
fuel_mass_fraction = [1]


# Film Cooling NASA
isnetropic = film.Isentropic(chamber_pressure, cea.T_static, cea.gamma)
mach = isnetropic.mach(geometry)[::-1]
T_aw_uncooled = isnetropic.adiabatic_wall_temp(mach, geometry, cea.Pr)
coolant = thermo.Mixture(['C2H5OH','H2O'], ws=[0.8,0.2], P=67e5, T=337)
film = film.FilmCooling(coolant, cea, 5.8, film_massflow, chamber_pressure, geometry[46,1], geometry)
T_aw_cooled = film.T_aw(film_start=46, film_end=90, mach=mach, T_aw_uncooled=T_aw_uncooled, n_holes=number_of_channels, chamber_pressure=chamber_pressure)[::-1]

heat = ht.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, ethanol90, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, cooling_channels, thermal_conductivity, method,T_aw_cooled)
heat.heatflux(geometry)

# Film cooling BETA
#coolant_velocity = 30
#film_start_idx = 48
#film = fm.FilmCooling(coolant, cea, total_massflow, chamber_pressure, coolant_velocity)
#film.film_effectiveness(film_start_idx, heat.section_length[::-1], heat.mach[::-1], heat.q[::-1], heat.q_rad[::-1], film_massflow, geometry)
#heat.heatflux(geometry, film.eta_arr[::-1], film_start_idx)

plt.rcParams.update({'font.size': 12})
f, axes = plt.subplots(4, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.coolant_pressure)
axes[3].set_ylabel('adiabatic wall temperature [K]')

plt.xlabel('x coordiante [m]')
plt.show()

print(max(heat.wall_temp))
print(max(heat.q)/1e6)
idx = np.where(heat.q == max(heat.q))
print('halpha gas: ', heat.halpha_gas[idx])
print('gas temp: ', heat.t_aw[idx])
print('halpha coolant: ', heat.halpha_coolant[idx])
print('coolant temp: ', heat.coolant_temp[idx])
print('radiation: ', heat.q_rad[idx])
'''
with open('temperature.csv', 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(["x_coordinate","y_coordinate","P_c","T_wall_inner","T_coolant","P_coolant"])
	for i in range(len(geometry[:,1])):
		writer.writerow([geometry[i,0], geometry[i,1], heat.P_chamber[len(geometry[:,1]) - i - 1], heat.wall_temp[len(geometry[:,1]) - i - 1], heat.coolant_temp[len(geometry[:,1]) - i - 1], heat.coolant_pressure[len(geometry[:,1]) - i - 1]])
'''

plt.plot(geometry[:,0][::-1], heat.q/1e6, label='total heat flux')
plt.plot(geometry[:,0][::-1], heat.halpha_gas*(heat.t_aw-heat.wall_temp)/1e6, label='convective heat flux')
plt.plot(geometry[:,0][::-1], heat.q_rad/1e6, label='radiation')
plt.grid()
plt.xlabel('x coordinate [m]')
plt.ylabel('heat flux [MW/m^2]')
plt.title('Pc 5 MPa, O/F 1.61, mu 0.18')
plt.legend(loc='best')
plt.show()
