import numpy as np
import thermo

import injectors as inj

geometry = np.genfromtxt('sparrow_08atm.txt', delimiter='', dtype=None, skip_header = 13) / 1000 
injection_idx = 46 
film_massflow = 0.411					# [kg/m^3]
coolant_temp = 337						# [K]
coolant_pressure = 6.7e6
pressure_drop = coolant_pressure - 5e6
gap_length = 4e-3
chamber_radius = geometry[injection_idx,1]


an_inj = inj.AnnulusInjector(['c2h5oh','h2o'], [0.8,0.2], coolant_temp, coolant_pressure, gap_length, 2*chamber_radius, film_massflow, pressure_drop)
an_inj.injector()
print(an_inj.mu)
print(an_inj.diameter*1000)

an_orifice = inj.AnnularOrifice(['c2h5oh','h2o'], [0.8,0.2], coolant_temp, coolant_pressure, gap_length, pressure_drop, film_massflow, chamber_radius)
an_orifice.injector()
print(an_orifice.mu)
print(an_orifice.diameter*1000)

n_holes = 42

liq_inj = inj.LiquidInjector(['c2h5oh','h2o'], [0.8,0.2], coolant_temp, coolant_pressure, gap_length, film_massflow/n_holes, pressure_drop, np.pi/2)
liq_inj.injector()
print(liq_inj.mu)
print(liq_inj.diameter*1000)
