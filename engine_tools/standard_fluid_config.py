import numpy as np
import engine_tools as et
import thermo
import rocketcea


#chamber parameters
chamber_pressure = 50e5 			# [Pa]
fuel_inlet_pressure = 75e5			# [Pa]
fuel_temperature = 288				# [K]
ox_temperature = 90 				# [K]
expansion_ratio = 7.93
pre_injection_pressure = 1.2*chamber_pressure	# [Pa]
fuel_injection_temperature = 453 	# [K]
chamber_diameter = 130e-3			# [m]
throat_diameter = 49.72e-3			# [m]

# CEA input values 
OF = 1.61		# 0.97 for injector face		
oxidiser = 'LOX'
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10])  # new fule blend for CEA
 

# Thermo input values 
total_massflow = 5.813 
fuel_massflow = total_massflow / (1+OF)
ox_massflow = total_massflow - fuel_massflow
fuel_composition = ['C2H5OH', 'H2O']
fuel_mass_fraction = [0.9, 0.1]
ox_composition = 'O2'



# Hot gas properties
# can choose beteen 'chamber'. 'throat' and 'exit' for cea.metric_cea_output
cea = et.CEA(ethanol90, oxidiser, chamber_pressure)
cea.metric_cea_output('chamber', OF, expansion_ratio)
#print(total_massflow*np.sqrt(cea.gamma*8314.5/cea.MW*cea.T_static))

# liquid properties
liquid_fuel = thermo.Mixture(fuel_composition, ws = fuel_mass_fraction, T=fuel_injection_temperature, P=pre_injection_pressure)
liquid_ox = thermo.Chemical(ox_composition, T=ox_temperature, P=pre_injection_pressure)

if __name__ == "__main__":
	print('gas static temperature: ', cea.T_static, '[K]')
	print('gas dynamic viscosity: ', cea.mu, '[Pa s]')
	print('gas Prandtl number: ', cea.Pr, '[-]')
	print('ratio of specific heats: ', cea.gamma, '[-]')
	print('gas thermal conductivity: ', cea.k, '[W/m/K]')
	#print('gas mole fractions: ', cea.mole_fractions)

	print('fuel density: ', liquid_fuel.rho, '[kg/m^3]')
	print('fuel dynamic visconsity', liquid_fuel.mu, '[Pa s]')
	print('fuel specific heat at constant pressure: ', liquid_fuel.Cpl, '[J/kg/K]')
	print('fuel thermal conductivity: ', liquid_fuel.k, '[W/m/K]')

	print('oxidiser density: ', liquid_ox.rho, '[kg/m^3]')
	print('oxidiser dynamic visconsity', liquid_ox.mu, '[Pa s]')
	print('oxidiser specific heat at constant pressure: ', liquid_ox.Cpl, '[J/kg/K]')
	print('oxidiser thermal conductivity: ', liquid_ox.k, '[W/m/K]')

