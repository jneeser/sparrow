import numpy as np
import thermo

import injectors
from standard_fluid_config import ox_massflow, fuel_massflow


fuel = thermo.Mixture(['c2h5oh', 'h2o'], ws=[0.8,0.2], T=288, P=7.5e6)
ox = thermo.Mixture(['o2'], [1], T=90, P=6.5e6)

# values form Wolfram Alpha
p_vap_water = 1204			#[Pa]
p_vap_ethanol = 3861 		#[Pa]
p_vap_ox = 97.18e3			#[Pa]

p_vap_fuel = fuel.zs[0] * p_vap_ethanol + fuel.zs[1] * p_vap_water			# molar averaged vapor pressure of fuel [Pa]

dp_fuel = fuel.P - p_vap_fuel
dp_ox = ox.P - p_vap_ox	


orifice_len = 3e-3
ethanol_bypass = injectors.LiquidInjector(['c2h5oh', 'h2o'], [0.8,0.2], fuel.T, fuel.P, orifice_len, fuel_massflow*0.1, dp_fuel, 0)
ethanol_bypass.injector()
lox_bypass = injectors.LiquidInjector(['o2'], [1], ox.T, ox.P, orifice_len, ox_massflow*0.1, dp_ox, 0)
lox_bypass.injector()

print(ethanol_bypass.mu)
print(ethanol_bypass.diameter*1000)

print(lox_bypass.mu)
print(lox_bypass.diameter*1000)