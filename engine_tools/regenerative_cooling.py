import numpy as np
import engine_tools as et
import thermo
from matplotlib import pyplot as plt

geometry = np.genfromtxt('sparrow_test_contour.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
wall_thickness = 1e-3
thermal_conductivity = 343
channel_hydrolic_diameter = np.ones(len(geometry[:,1]))*1e-3

chamber_pressure = 40e5 	
coolant_inlet_pressure = 60e5
coolant_temperature = 288
OF = 1.4472			# reduced OF ratio to reflect 90% C2H5OH
total_massflow = 5.5 
coolant = ['H2O']
coolant_mass_fraction = [1]
coolant_massflow = 2.2
fuel = 'C2H5OH'
oxidiser = 'LOX'


heat = et.Heattransfer(coolant, coolant_mass_fraction, coolant_massflow, total_massflow, fuel, oxidiser, OF, chamber_pressure, coolant_temperature, coolant_inlet_pressure, geometry, wall_thickness, thermal_conductivity)
#heat = et.Heattransfer(['CH4'], [1], 1.32, 5.5, 'CH4', 'LOX', 3.16, 40e5, 110, 60e5, geometry, wall_thickness, thermal_conductivity)

heat.heatflux(channel_hydrolic_diameter, geometry, 1000, False)

f, axes = plt.subplots(5, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.coolant_pressure1)
axes[3].set_ylabel('coolant pressure [MPa]')

axes[4].plot(geometry[:,0][::-1], heat.adiabatic_wall_temperature)
axes[4].set_ylabel('coolant temperature [K]')

plt.xlabel('x coordiante [m]')
plt.show()

print(max(heat.wall_temp))
