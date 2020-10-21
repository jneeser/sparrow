from matplotlib import pyplot as plt
import numpy as np


data = np.genfromtxt('optimised_geometry.csv', delimiter=',', dtype=None, skip_header = 1)
geometry = np.genfromtxt('sparrow_contour_1_5.txt', delimiter='', dtype=None, skip_header = 13) / 1000 					# conversion to [m]

x = geometry[:,0] 
y = geometry[:,1]
halpha_gas_arr = data[:,2]
q_total_arr = data[:,3]
wall_t_arr = data[:,4]
Rei_arr = data[:,5]
Reo_arr = data[:,6]
dp_arr = data[:,7]
sec_length_arr = data[:,8]
dhi_arr = data[:,9]
dho_arr = data[:,10]
wt1_arr = data[:,11]
wto_arr = data[:,12]
rf1i_arr = data[:,13]
rf1o_arr = data[:,14]
rf2_arr = data[:,15]
t_arr = data[:,16]
hi_arr = data[:,17]
ho_arr = data[:,18]
stress_ratio_arr = data[:,19]

f, axes = plt.subplots(5, 1)
axes[0].plot(geometry[:,0], geometry[:,1]*1000)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(geometry[:,0], wall_t_arr)
axes[1].set_ylabel('wall temperature [K]')

axes[2].plot(geometry[:,0], dp_arr/1e6)
axes[2].set_ylabel('Pressure drop [MPa/m]')

axes[3].plot(geometry[:,0], stress_ratio_arr)
axes[3].set_ylabel('stress ratio [-]')

axes[4].plot(geometry[:,0], q_total_arr/1e6)
axes[4].set_ylabel('heat flux [Mw/m^2]')

plt.xlabel('x coordiante [m]')
plt.show()

print('total pressure drop: ', sum(dp_arr*sec_length_arr/1e6), '[MPa]')