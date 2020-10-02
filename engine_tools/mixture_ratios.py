import numpy as np
import engine_tools as et
import rocketcea
from rocketcea.cea_obj import CEA_Obj
import pylab as pl


chamber_pressure = 50e5 			# [Pa]
expansion_ratio = 8


# CEA input values 
oxidiser = 'LOX'
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])  # new fule blend for CEA
 


from rocketcea.cea_obj import CEA_Obj

pcM = [40e5, 50e5, 60e5]

ispObj = CEA_Obj(propName='', oxName='LOX', fuelName=ethanol90)

for Pc in pcM:
	ispArr = []
	MR = 1.2
	mrArr = []
	while MR < 2:
		ispArr.append(ispObj.get_Isp(Pc=Pc*0.000145038, MR=MR, eps=expansion_ratio))
		mrArr.append(MR)
		MR += 0.001
	pl.plot(mrArr, ispArr, label='Pc=%g Pa'%Pc)

pl.legend(loc='best')
pl.grid(True)
pl.title( ispObj.desc )
pl.xlabel( 'Mixture Ratio [-]' )
pl.ylabel( 'Isp [s]' )
pl.show()


idx = np.where(np.array(ispArr) == max(ispArr))[0][0]
print('maximum Isp mixture ratio: ', mrArr[idx])
