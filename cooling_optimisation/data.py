from matplotlib import pyplot as plt
import numpy as np
import csv

data = np.genfromtxt('optimised_geometry.csv', delimiter=',', dtype=None, skip_header = 1)

x = data[:,0] 
y = data[:,1]
halpha_gas_arr = data[:,2]
q_total_arr = data[:,3]
wall_t_arr = data[:,4]
tbc_t_arr = data[:,5]
Rei_arr = data[:,6]
Reo_arr = data[:,7]
dp_arr = data[:,8]
sec_length_arr = data[:,9]
dhi_arr = data[:,10]
dho_arr = data[:,11]
wt1_arr = data[:,12]
wto_arr = data[:,13]
rf1i_arr = data[:,14]
rf1o_arr = data[:,15]
rf2_arr = data[:,16]
t_arr = data[:,17]
hi_arr = data[:,18]
ho_arr = data[:,19]
stress_ratio_arr = data[:,20]

f, axes = plt.subplots(5, 1)
axes[0].plot(x, y*1e3)
axes[0].set_ylabel('contour height [mm]')

axes[1].plot(x, t_arr*1e3)
axes[1].set_ylabel('t [mm]')

axes[2].plot(x, np.round(wto_arr*1e3,4))
axes[2].set_ylabel('Wto [mm]')

axes[3].plot(x, np.round(rf2_arr*1e3,4))
axes[3].set_ylabel('R2 [mm]')

axes[4].plot(x, np.round(rf1o_arr*1e3,4))
axes[4].set_ylabel('R1o [mm]')

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

