from matplotlib import pyplot as plt
import numpy as np
import csv
from scipy import signal

data = np.genfromtxt('optimised_geometry.csv', delimiter=',', dtype=None, skip_header = 1)

x = data[:,0] 
y = data[:,1]
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

b, a = signal.butter(1, 0.1)
zi = signal.lfilter_zi(b, a)
z, _ = signal.lfilter(b, a, rf1i_arr, zi=zi*rf1i_arr[0])

f, axes = plt.subplots(5, 1)
axes[0].plot(x, y*1e-3)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(x, wall_t_arr)
axes[1].set_ylabel('wt1 [mm]')

axes[2].plot(x, rf2_arr*1e3)
axes[2].set_ylabel('rf2 [mm]')

axes[3].plot(x, rf1i_arr*1e3)
axes[3].set_ylabel('rf1i [mm]')

axes[4].plot(x, z*1e3)
axes[4].set_ylabel('rf1o filtered [mm]')

plt.xlabel('x coordinate [m]')
plt.show()

print('total pressure drop: ', sum(dp_arr*sec_length_arr/1e6), '[MPa]')
print(wt1_arr[0])


# export geometric parameters in mm for catia import 
with open('geometrical_params.csv', 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(["wt1", "wto", "rf1i", "rf1o", "rf2", "t"]
	)
	for i in range(len(y)):
		if i%2 == 0:
			idx = i
			writer.writerow([np.round(wt1_arr[idx]*1e3,3), np.round(wto_arr[idx]*1e3,3), np.round(rf1i_arr[idx]*1e3,3), np.round(rf1o_arr[idx]*1e3,3), np.round(rf2_arr[idx]*1e3,3), np.round(t_arr[idx]*1e3,3)
			])

