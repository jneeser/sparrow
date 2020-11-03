import numpy as np
import engine_tools as et
import thermo
import rocketcea
from matplotlib import pyplot as plt
import csv

data = np.genfromtxt('optimised_geometry.csv', delimiter=',', dtype=None, skip_header = 1)
geometry = np.genfromtxt('sparrow_contour_1_5.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
wall_thickness = data[:,11]
wall_thickness_tbc = 0.1e-3
thermal_conductivity = 24							# Inconel at 800 C
thermal_conductivity_tbc = 0.8
channel_hydrolic_diameter = (data[:,9] + data[:,10])/2
number_of_channels = 42*2

# global properties
chamber_pressure = 50e5 			# [Pa]
fuel_inlet_pressure = 70e5			# [Pa]
fuel_temperature = 288				# [K]
ox_temperature = 90 				# [K]
expansion_ratio = 7.93

# CEA input values 
OF = 1.61					# actual OF ratio 
oxidiser = 'LOX'
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])  # new fuel blend for CEA
 

# Thermo input values 
total_massflow = 5.813 						# [kg/s]
fuel_massflow = total_massflow / (1+OF)
ox_massflow = total_massflow - fuel_massflow
fuel_composition = ['C2H5OH']
fuel_mass_fraction = [1]


heat = et.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, ethanol90, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, number_of_channels, thermal_conductivity, thermal_conductivity_tbc, wall_thickness_tbc)
#heat = et.Heattransfer(['CH4'], [1], 1.32, 5.5, 'CH4', 'LOX', 3.16, 40e5, 110, 60e5, geometry, wall_thickness, thermal_conductivity)

heat.heatflux(channel_hydrolic_diameter, geometry, wall_thickness, 1250, False)

f, axes = plt.subplots(5, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.coolant_pressure)
axes[3].set_ylabel('coolant pressure [MPa]')

axes[4].plot(geometry[:,0][::-1], heat.coolant_temp)
axes[4].set_ylabel('coolant temperature [K]')

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

