import numpy as np
import engine_tools as et
import thermo
from matplotlib import pyplot as plt

geometry = np.genfromtxt('sparrow_test_contour.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
wall_thickness = 1e-3
thermal_conductivity = 343
channel_hydrolic_diameter = np.ones(len(geometry[:,1]))*1e-3

chamber_pressure = 40e5 				# [Pa]
coolant_inlet_pressure = 60e5			# [Pa]
coolant_temperature = 288				# [K]
OF = 1.6
total_massflow = 5.5 					# [kg/s]
t_burn = 40								# [s]
density_graphite = 1700					# [kg/m^3]	

heat = et.Heattransfer(['H2O'], [1], 2.2, total_massflow, 'C2H5OH', 'LOX', OF, chamber_pressure, coolant_temperature, coolant_inlet_pressure, geometry, wall_thickness, thermal_conductivity)

y = geometry[:,1][::-1]
x = geometry[:,0][::-1]
ablator_thickness = np.ndarray(len(y))
local_ablator_mass = np.ndarray(len(y))

initial_guess = np.ndarray(len(y))					
diverging = True
for i in range(len(y)):
	initial_guess[i] = diverging
	if abs(y[i] - heat.throat_diameter/2) < 1e-6:
		diverging = False

local = et.Isentropic(heat.chamber_pressure, heat.cea.T_static, heat.cea.gamma)

# Iterate over each chamber lcoation 
for i in range(len(y)):
	if i == 0:
		section_length = 0
	else:
	    section_length = np.sqrt((x[i] - x[i-1])**2 + (y[i]-y[i-1])**2)

	local_area = np.pi*y[i]**2
	mach = local.mach(local_area, np.pi*heat.throat_diameter**2/4, initial_guess[i])
	adiabatic_wall_temperature = heat.adiabatic_wall_temp(mach, initial_guess[i])

	c = 1.05
	Rr = 0.3
	Rv = 0.41
	cp = 0.38			# [Btu/lb/F]
	Lp = 686			# [Btu/lb]
	Td = 1460			# [R]
	pc = 580			# [psi]
	rho = 0.061			# [lb/in^3]
	k = 9.8e-6			# [Btu/in^2/F/in]
	t_burn = 40			# [s]
	taw_imperial = adiabatic_wall_temperature*9/5

	a = c*(2*k*t_burn / (Rr*Rv*cp*rho) * np.log(1 + Rr*Rv*cp*(taw_imperial - Td)/Lp))**0.5 * (pc/100)**0.4
	ablator_thickness[i] = a*0.0254
	local_ablator_mass[i] = np.pi*y[i]*2*ablator_thickness[i]

print('total ablator mass: ', sum(local_ablator_mass) , '[kg]')

f, axes = plt.subplots(2, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], ablator_thickness*1000)
axes[1].set_ylabel('ablator thickness [mm]')

plt.xlabel('x coordiante [m]')
plt.show()
