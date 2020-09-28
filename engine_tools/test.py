import numpy as np
import engine_tools_channels as ft
import thermo
from matplotlib import pyplot as plt

geometry = np.genfromtxt('sparrow_test_contour.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
wall_thickness = 1e-3
thermal_conductivity = 343
channel_hydrolic_diameter = np.ones(len(geometry[:,1]))*1.5e-3
channel_width = 2e-3
channel_height = 1e-3
n_channels = 72


heat = ft.Heattransfer(['C2H5OH'], [1], 2.2, 5.5, 'C2H5OH', 'LOX', 1.6, 40e5, 288, 60e5, geometry, wall_thickness, thermal_conductivity)
#heat = ft.Heattransfer(['CH4'], [1], 1.32, 5.5, 'CH4', 'LOX', 3.16, 40e5, 110, 60e5, geometry, wall_thickness, thermal_conductivity)

heat.heatflux(channel_width, channel_height, n_channels, geometry, 1000, True)

f, axes = plt.subplots(5, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], heat.wall_temp)
axes[1].set_ylabel('Wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], heat.q/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], heat.coolant_pressure/1e6)
axes[3].set_ylabel('coolant pressure [MPa]')

axes[4].plot(geometry[:,0][::-1], heat.adiabatic_wall_temperature)
axes[4].set_ylabel('coolant temperature [K]')
plt.xlabel('x coordiante [m]')
plt.show()



print(max(heat.wall_temp))
