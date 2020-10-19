import heatflux
import standard_fluid_config as std

import numpy as np
import thermo
import rocketcea
from matplotlib import pyplot as plt
import csv


geometry = np.genfromtxt('sparrow_contour_1_5.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
isnetropic = heatflux.Isentropic(std.chamber_pressure, std.cea.T_static, std.cea.gamma)
mach = isnetropic.mach(geometry)
t_aw = isnetropic.adiabatic_wall_temp(mach, geometry, std.cea.Pr)


wall_thickness = np.ones(len(geometry[:,1]))*0.4e-3
thermal_conductivity = 28							# Inconel at 800 C
channel_hydrolic_diameter = np.ones(len(geometry[:,1]))*1.4e-3
number_of_channels = 60


# coolant properties 
fuel_composition = ['C2H5OH']
fuel_mass_fraction = [1]


heat = heatflux.Heattransfer(fuel_composition, fuel_mass_fraction, std.fuel_massflow, std.total_massflow, std.cea, std.chamber_pressure, std.fuel_temperature, std.fuel_inlet_pressure, geometry, number_of_channels, thermal_conductivity)

heat.heatflux(channel_hydrolic_diameter, geometry, wall_thickness, mach, t_aw)

f, axes = plt.subplots(4, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.coolant_temp)
axes[3].set_ylabel('coolant temperature [K]')

plt.xlabel('x coordiante [m]')
plt.show()

print(max(heat.wall_temp))