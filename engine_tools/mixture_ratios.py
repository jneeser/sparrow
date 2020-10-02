import numpy as np
import engine_tools as et
import thermo
import rocketcea
from matplotlib import pyplot as plt


chamber_pressure = 50e5 			# [Pa]
expansion_ratio = 7.93


# CEA input values 
oxidiser = 'LOX'
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])  # new fule blend for CEA
 

mixture_ratios = np.arange(1, 2, 0.001)
cstar = np.ndarray(len(mixture_ratios))
t_static = np.ndarray(len(mixture_ratios))
isp = np.ndarray(len(mixture_ratios))

for i in range(len(mixture_ratios)):
	cea = et.CEA(ethanol90, oxidiser, chamber_pressure)
	cea.metric_cea_output('chamber', mixture_ratios[i], expansion_ratio)
	cstar[i] = cea.cstar
	t_static[i] = cea.T_static
	isp[i] = cea.isp


idx = np.where(cstar == max(cstar))
print('maximum cstar mixture ratio: ', mixture_ratios[idx])


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Mixture Ratio [-]')
ax1.set_ylabel('cstar [m/s]', color=color)
ax1.plot(mixture_ratios, cstar, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Static Temperature [K]', color=color)  # we already handled the x-label with ax1
ax2.plot(mixture_ratios, t_static ,color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.show()


