import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

import injectors as inj

def filtering(data, order, cutoff):
	# lowpass butterworth filter applied three times
	b, a = signal.butter(order, cutoff)

	zi = signal.lfilter_zi(b, a)
	z1,_ = signal.lfilter(b, a, data, zi=zi*data[0])
	z2,_ = signal.lfilter(b, a, z1, zi=zi*z1[0])
	z3,_ = signal.lfilter(b, a, z2, zi=zi*z2[0])

	return z3

def get_ps_data(col, dataset):
	slope = 6000
	offset = -23
	return dataset[:,col]*slope + offset 

def get_lc_data(col,dataset):
	slope = -25070
	offset = -21.27
	return dataset[:,col]*slope + offset 

def get_time(col, dataset, frequency):
	slope = 1/frequency
	offset = dataset[:,col][0]
	return (dataset[:,col] - offset)*slope/1000 

def get_mass_flow(mass, time, intervall):
	def func(x, a, b):
		return a*x +b

	mass_flows = []
	times = []
	idx = 0 

	while idx < len(mass):
		popt,_ = curve_fit(func, time[idx:idx+intervall], mass[idx:idx+intervall])
		mass_flows.append(popt[0])
		times.append(time[idx+int(intervall/2)])
		idx += intervall

	return np.array(mass_flows), np.array(times)




raw_ps_data = np.genfromtxt('Expulsion_Test_1_Pressure_Data.csv', delimiter=',', skip_header=77500, max_rows=5000)
raw_lc_data = np.genfromtxt('Expulsion_Test_1_Loadcell_Data.csv', delimiter=',', skip_header=77500, max_rows=5000) 

ps_time = get_time(8, raw_ps_data, 1000)
ps_tank = get_ps_data(0,raw_ps_data)
ps_tank_filt = filtering(ps_tank, 3, 0.04)
ps_manifold = get_ps_data(1,raw_ps_data)
ps_manifold_filt = filtering(ps_manifold, 3, 0.04)


lc_time = get_time(4, raw_lc_data, 1000)
lc_tank = get_lc_data(1,raw_lc_data)
lc_tank_filt = filtering(lc_tank, 3, 0.02)


mass_flow, mass_flow_times = get_mass_flow(lc_tank, lc_time, 1000)
#mass_flow = filtering(mass_flow, 3, 0.08)


plt.plot(ps_time, ps_manifold_filt, label='manifold pressure')
plt.plot(ps_time, ps_tank_filt, label='tank pressure')
plt.legend(loc='best')
plt.show()

plt.plot(lc_time, lc_tank_filt)
plt.show()

plt.scatter(mass_flow_times, mass_flow)
plt.show()