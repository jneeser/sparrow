import heatflux
import standard_fluid_config as std
import pressuredrop 

import numpy as np
import thermo
import rocketcea
from matplotlib import pyplot as plt
import csv
import scipy.optimize
from scipy import interpolate
import time


# coolant properties 
fuel_composition = ['C2H5OH']
fuel_mass_fraction = [1]
number_of_channels = 42*2
wall_thickness_tbc = 0.1e-3
thermal_conductivity = 28							# Inconel at 800 C
thermal_conductivity_tbc = 2

# create initial objects 
geometry = np.genfromtxt('sparrow_contour_coarse.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]
isnetropic = heatflux.Isentropic(std.chamber_pressure, std.cea.T_static, std.cea.gamma)
mach = isnetropic.mach(geometry)
t_aw = isnetropic.adiabatic_wall_temp(mach, geometry, std.cea.Pr)

heat = heatflux.Heattransfer(fuel_composition, fuel_mass_fraction, std.fuel_massflow, std.total_massflow, std.cea, std.chamber_pressure, std.fuel_temperature, std.fuel_inlet_pressure, geometry, number_of_channels, thermal_conductivity, thermal_conductivity_tbc, wall_thickness_tbc)

y_coordinates = geometry[:,1][::-1]
x_coordinates = geometry[:,0][::-1]

###############################
# starting the section iterator 
###############################

tol = 1e-2 							# acceptable temeprature difference 
max_iter = 1000

# empty output arrays for pressuredrop.py output
fi_arr = np.ndarray(len(y_coordinates))
dp_arr = np.ndarray(len(y_coordinates))
vi_arr = np.ndarray(len(y_coordinates))
Rei_arr = np.ndarray(len(y_coordinates))
dhi_arr = np.ndarray(len(y_coordinates))
dho_arr = np.ndarray(len(y_coordinates))
vo_arr = np.ndarray(len(y_coordinates))
Reo_arr = np.ndarray(len(y_coordinates))
hi_arr = np.ndarray(len(y_coordinates))
ho_arr = np.ndarray(len(y_coordinates))
wt1_arr = np.ndarray(len(y_coordinates))
rf1i_arr = np.ndarray(len(y_coordinates))
rf1o_arr = np.ndarray(len(y_coordinates))
rf2_arr = np.ndarray(len(y_coordinates))
wto_arr = np.ndarray(len(y_coordinates))
t_arr = np.ndarray(len(y_coordinates))

# empty arrays for heatflux outputs 
wall_t_arr = np.ndarray(len(y_coordinates))
coolant_t_arr = np.ndarray(len(y_coordinates))
coolant_p_arr = np.ndarray(len(y_coordinates))
q_total_arr = np.ndarray(len(y_coordinates))
q_rad_arr = np.ndarray(len(y_coordinates))
halpha_gas_arr = np.ndarray(len(y_coordinates))
sec_length_arr = np.ndarray(len(y_coordinates))

# material 
in718 = pressuredrop.metal(E=200e9,k=20,v=0.33,alpha=12e-6,sig_yield=([1150e6,1150e6,950e6,650e6],[273,50+273,700+273,850+273]))


#TODO implement pressure drops 
# initial guess
wall_temperature = 600
coolant_temp = 400
radiation = 0
coolant_pressure = 70e5
halpha = 3000
t = 1.2e-3
wt1 = 0.6e-3
wt2 = 0.4e-3
rf1 = 0.1e-3
rf2 = 0.1e-3

t1 = time.time()
# optimisation 
for i in range(len(y_coordinates)):

	if i == 0:
		section_length = 0
	else:
		section_length = np.sqrt((x_coordinates[i] - x_coordinates[i-1])**2 + (y_coordinates[i]-y_coordinates[i-1])**2)
	
	# iteration parameters 
	iteration = 0
	difference = 1

	while difference > tol:
		initial_params = pressuredrop.parameters(ri=y_coordinates[i],t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=int(number_of_channels/2))
		initial_heat = pressuredrop.sim(t_aw[i], heat.coolant.T,halpha,radiation,coolant_pressure,std.chamber_pressure,0)
		try:
			new_params = pressuredrop.physics(initial_params,in718,initial_heat,heat.coolant)
		except:
			print('skipped section: ', i, ' due to infeasible solution')
			pass
		avg_hydrolic_diameter = (new_params.dhi + new_params.dho) / 2 
		wall_thickness = new_params.par.wt1

		heat_flux, new_wall_temp, tbc_wall_temp, Re, Nu, radiation, halpha = heat.iterator(y_coordinates[i], avg_hydrolic_diameter, section_length, wall_thickness, mach[i], t_aw[i])

		difference = abs(new_wall_temp - wall_temperature)

		iteration += 1
		print('section: ', i, ' sub-iteration: ', iteration, ':')
		print('	temperature difference: ', difference)
		print('	wall temperature: ', wall_temperature)
		print('	coolant temperature: ', heat.coolant.T)

		if iteration > max_iter:
			raise ValueError('Non-convergence, iteration number exceeded ', max_iter)

		wall_temperature = new_wall_temp	
		wt1 = new_params.par.wt1
		rf1 = new_params.par.rf1i
		rf2 = new_params.par.rf2
		t = new_params.par.t
		wt2 = new_params.par.wt2

	T_new = heat.coolant.T + heat_flux*2*np.pi*y_coordinates[i]*section_length / (heat.coolant_massflow*heat.coolant.Cp) 
	heat.coolant.calculate(P=heat.coolant.P, T=T_new)


	# updating arrays 
	fi_arr[i] = new_params.fi
	dp_arr[i] = new_params.dp
	vi_arr[i] = new_params.vi
	Rei_arr[i] = new_params.Rei
	dhi_arr[i] = new_params.dhi
	dho_arr[i] = new_params.dho
	vo_arr[i] = new_params.vo
	Reo_arr[i] = new_params.Reo
	hi_arr[i] = new_params.hi
	ho_arr[i] = new_params.ho
	wt1_arr[i] = new_params.par.wt1
	rf1i_arr[i] = new_params.par.rf1i
	rf1o_arr[i] = new_params.par.rf1o
	rf2_arr[i] = new_params.par.rf2
	wto_arr[i] = new_params.wto
	t_arr[i] = new_params.par.t

	# empty arrays for heatflux outputs 
	wall_t_arr[i] = wall_temperature
	coolant_t_arr[i] = T_new
	q_total_arr[i] = heat_flux
	q_rad_arr[i] = radiation
	halpha_gas_arr[i] = halpha
	sec_length_arr[i] = section_length

t2 = time.time()
print('optimisation runtime: ', t2-t1, '[s]')

f, axes = plt.subplots(4, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0][::-1], wall_t_arr)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0][::-1], q_total_arr/1e6)
axes[2].set_ylabel('heat flux [MW/m^2]')

axes[3].plot(geometry[:,0][::-1], dp_arr)
axes[3].set_ylabel('coolant pressure [MPa]')

plt.xlabel('x coordiante [m]')
plt.show()


with open('optimised_geometry.csv', 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(["x coordinate","y_coordinate","gas heat transfer coefficient", "heat flux", "max wall temperature", "Re inner", "Re outer", "pressure drops", "section lenghts","hydrolic diameter inner",
					 "hydrolic diameter outer", "wt1", "wto", "rf1 inner", "rf1 outer", "rf2", "inner wall thickness", "hi", "ho"]
	)
	for i in range(len(geometry[:,1])):
		idx = len(geometry[:,1]) - i - 1
		writer.writerow([geometry[i,0], geometry[i,1], halpha_gas_arr[idx], q_total_arr[idx], wall_t_arr[idx], Rei_arr[idx], Reo_arr[idx], dp_arr[idx], sec_length_arr[idx], dhi_arr[idx],
						dho_arr[idx], wt1_arr[idx], wto_arr[idx], rf1i_arr[idx], rf1o_arr[idx], rf2_arr[idx], t_arr[idx], hi_arr[idx], ho_arr[idx]
		])

