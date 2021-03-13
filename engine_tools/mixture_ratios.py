import numpy as np
import engine_tools as et
import rocketcea
from rocketcea.cea_obj import CEA_Obj
import pylab as pl
from matplotlib import pyplot as plt

import standard_fluid_config as std


pcM = [50e5]
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[80,20])  # new fule blend for CEA
ispObj = CEA_Obj(propName='', oxName='LOX', fuelName=ethanol90)

for Pc in pcM:
	ispArr = []
	TArr = []
	MR = 0.6
	mrArr = []
	while MR < 2:
		ispArr.append(ispObj.get_Isp(Pc=Pc*0.000145038, MR=MR, eps=std.expansion_ratio))
		TArr.append(ispObj.get_Tcomb(Pc=Pc*0.000145038, MR=MR)*0.555556)
		mrArr.append(MR)
		MR += 0.001
	pl.plot(mrArr, ispArr, label='Pc=%g Pa'%Pc)

pl.legend(loc='best')
pl.grid(True)
pl.title( ispObj.desc )
pl.xlabel( 'Mixture Ratio [-]' )
pl.ylabel( 'Isp [s]' )
pl.show()

pl.plot(mrArr, TArr, label='Pc=%g Pa'%Pc)
pl.legend(loc='best')
pl.grid(True)
pl.title( ispObj.desc )
pl.xlabel( 'Mixture Ratio [-]' )
pl.ylabel( 'T [K]' )
pl.show()

    
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Mixture Ratio [-]')
ax1.set_ylabel('Isp [s]', color=color)
ax1.plot(mrArr, ispArr, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Temperature [K]', color=color)  # we already handled the x-label with ax1
ax2.plot(mrArr, TArr, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.show()




idx = np.where(np.array(ispArr) == max(ispArr))[0][0]
print('maximum Isp mixture ratio: ', mrArr[idx])
