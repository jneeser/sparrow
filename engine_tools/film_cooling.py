import numpy as np
import thermo
import rocketcea
from matplotlib import pyplot as plt
import csv

import heat_transfer as ht

data = np.genfromtxt('optimised_geometry_50bar.csv', delimiter=',', dtype=None, skip_header = 1)
geometry = np.genfromtxt('sparrow_50bar.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
wall_thickness = data[:,12][::-1]
wall_thickness_tbc = 0
thermal_conductivity = 24							# Inconel at 800 C
thermal_conductivity_tbc = 1.2
channel_hydraulic_diameter = data[:,10][::-1]
number_of_channels = 42*2
method = 'cinjarew'

# global properties
chamber_pressure = 50e5 			# [Pa]
fuel_inlet_pressure = 75e5			# [Pa]
fuel_temperature = 288				# [K]
ox_temperature = 90 				# [K]
expansion_ratio = 12

# CEA input values 
OF = 1.61	
oxidiser = 'LOX'
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])  # new fuel blend for CEA
 

# Thermo input values 
total_massflow = 5.8 						# [kg/s]
fuel_massflow = total_massflow / (1+OF)
ox_massflow = total_massflow - fuel_massflow
fuel_composition = ['C2H5OH']
fuel_mass_fraction = [1]


heat = ht.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, ethanol90, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, number_of_channels, thermal_conductivity, method, thermal_conductivity_tbc, wall_thickness_tbc)

heat.heatflux(channel_hydraulic_diameter, geometry, wall_thickness, 1250, False)

plt.rcParams.update({'font.size': 12})
f, axes = plt.subplots(4, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.tbc_wall_temp)
axes[3].set_ylabel('coolant temperature [K]')

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
wt1_arr = data[:,11][::-1]
wto_arr = data[:,12][::-1]
rf1i_arr = data[:,13][::-1]
rf1o_arr = data[:,14][::-1]
rf2_arr = data[:,15][::-1]
t_arr = data[:,16][::-1]

idx = np.where(heat.q == max(heat.q))
print(heat.coolant_temp[idx])
print(heat.flowvelocity[idx])
print(wt1_arr[idx])
print(wto_arr[idx])
print(rf1i_arr[idx])
print(rf1o_arr[idx])
print(rf2_arr[idx])
print(t_arr[idx])


plt.plot(geometry[:,0][::-1], heat.q/1e6, label='total heat flux')
plt.plot(geometry[:,0][::-1], heat.halpha_gas*(heat.t_aw-heat.tbc_wall_temp)/1e6, label='convective heat flux')
plt.plot(geometry[:,0][::-1], heat.q_rad/1e6, label='radiation')
plt.grid()
plt.xlabel('x coordinate [m]')
plt.ylabel('heat flux [MW/m^2]')
plt.legend(loc='best')
plt.show()
