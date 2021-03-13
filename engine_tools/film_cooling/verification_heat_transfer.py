import numpy as np
import thermo
import rocketcea
from matplotlib import pyplot as plt
from scipy import interpolate

import heat_transfer as ht
import geom_class as gc

geometry = np.genfromtxt('rl-10.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
t = 2.2e-3
wt1 = 0.33e-3
wt2 = 0.33e-3
rf1 = 0.2e-3
rf2 = 0.2e-3
number_of_channels = 90
cooling_channels = gc.parameters(ri=130e-3,t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=number_of_channels)
#cooling_channels.cooling_geometry(geometry[:,1][::-1])
cooling_channels.dhi_arr = np.ones(len(geometry[:,1]))*2.6e-3
cooling_channels.wt1_arr = np.ones(len(geometry[:,1]))*0.33e-3
cooling_channels.Ai_arr = cooling_channels.dhi_arr**2*np.pi/4


thermal_conductivity = 386							# copper at 600 C
method = 'cinjarew'

# global properties
chamber_pressure = 33.09e5 			# [Pa]
fuel_inlet_pressure = 130e5			# [Pa]
fuel_temperature = 64				# [K]
ox_temperature = 90 				# [K]
expansion_ratio = 61

# CEA input values 
OF = 5.26	
oxidiser = 'LOX'
fuel = 'H2'

# Thermo input values 
total_massflow = 17.06						# [kg/s]
fuel_massflow = 2.709
fuel_composition = ['H2']
fuel_mass_fraction = [1]

heat = ht.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, fuel, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, cooling_channels, thermal_conductivity, method)

heat.heatflux(geometry)

plt.rcParams.update({'font.size': 12})
f, axes = plt.subplots(4, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.coolant_pressure)
axes[3].set_ylabel('coolant temperature [K]')

plt.xlabel('x coordiante [m]')
plt.show()

print(max(heat.wall_temp))

method2 = 'standard-bartz'
method3 = 'modified-bartz'
heat2 = ht.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, fuel, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, cooling_channels, thermal_conductivity, method2)
heat3 = ht.Heattransfer(fuel_composition, fuel_mass_fraction, fuel_massflow, total_massflow, fuel, oxidiser, OF, chamber_pressure, fuel_temperature, fuel_inlet_pressure, geometry, cooling_channels, thermal_conductivity, method3)
heat2.heatflux(geometry)
heat3.heatflux(geometry)


throat_idx = np.where(geometry[:,1] == min(geometry[:,1]))[0][0]
q_exp_arr = np.array([6.2,6.3,6.2,7,10,15.8,19.9,16.7,10.2,5.8,3.7,2.2,1.7,1.8,1,0.8,0.7])*1.64 								# conversion to [MW/m^2]
pos_arr = np.array([-12,-10,-8,-6,-4,-2,0,1,2.7,4.5,7,9,16,17,25,33,42])*0.0254 + np.ones(len(q_exp_arr))*geometry[throat_idx,0]	# correcting position and conversion to [m]


plt.plot(geometry[:,0][::-1], heat.q/1e6, label='heat flux cinjarew')
plt.plot(geometry[:,0][::-1], heat2.q/1e6, label='heat flux SB', linestyle='dashed')
plt.plot(geometry[:,0][::-1], heat3.q/1e6, label='heat flux MB', linestyle='dotted')
plt.scatter(pos_arr, q_exp_arr, color='red', label='reference heat flux')
plt.grid()
plt.xlabel('x coordinate [m]')
plt.ylabel('heat flux [MW/m^2]')
plt.title('RL-10 heat flux')
plt.legend(loc='best')
plt.show()

q = interpolate.interp1d(heat.q,geometry[:,0][::-1],fill_value='extrapolate')
error = []
for i in range(len(q_exp_arr)):
	error.append((q_exp_arr[i]- q(pos_arr[i]))/1e6 / q_exp_arr[i])

print(abs(np.average(np.array(error))))