import numpy as np

from analysis import filtering, get_ps_data, get_lc_data, get_mass_flow, get_discharge_coeff, get_time, get_dyn_pressure
from plotting import plot_spectrum, pressure_plot, loadcell_plot, massflow_plot


# Data import 
raw_ps_data = np.genfromtxt('Hot2PS.csv', delimiter=',', skip_header=1)
raw_lc_data = np.genfromtxt('Hot2LC.csv', delimiter=',', skip_header=1) 

print(len(raw_ps_data[:,6]) / abs(raw_ps_data[:,6][0] - raw_ps_data[:,6][-1])/1e6)

# Filtering parameters
ps_filter_order = 5
ps_filter_cutoff = 0.08
lc_filter_order = 3
lc_filter_cutoff = 0.06

# Data processing 
ps_time = get_time(6, raw_ps_data, 1e6, -3065)
ps_ethanol_tank = get_ps_data(2,raw_ps_data, 6250, -24)
ps_ethanol_tank_filt = filtering(ps_ethanol_tank, ps_filter_order, ps_filter_cutoff)

ps_lox_tank = get_ps_data(0,raw_ps_data, 6250, -24)
ps_lox_tank_filt = filtering(ps_lox_tank, ps_filter_order, ps_filter_cutoff)

ps_ethanol_manifold = get_ps_data(3,raw_ps_data, 6250, -24)
ps_ethanol_manifold_filt = filtering(ps_ethanol_manifold, ps_filter_order, ps_filter_cutoff)

ps_lox_manifold = get_ps_data(1,raw_ps_data, 15100, -59)
ps_lox_manifold_filt = filtering(ps_lox_manifold, ps_filter_order, ps_filter_cutoff)

ps_chamber_1 = get_ps_data(4,raw_ps_data, 3750, -14)
ps_chamber_1_filt = filtering(ps_chamber_1, ps_filter_order, ps_filter_cutoff)

ps_chamber_2 = get_ps_data(5,raw_ps_data, 3750, -14)
ps_chamber_2_filt = filtering(ps_chamber_2, ps_filter_order, ps_filter_cutoff)


lc_time = get_time(3, raw_lc_data, 1e6, -3065)
lc_ethanol_tank = get_lc_data(2,raw_lc_data, 169520, -156)
lc_ethanol_tank_filt = filtering(lc_ethanol_tank, lc_filter_order, lc_filter_cutoff)

lc_lox_tank = get_lc_data(0,raw_lc_data, 453283, -164)
lc_lox_tank_filt = filtering(lc_lox_tank, lc_filter_order, lc_filter_cutoff)

lc_thrust = get_lc_data(1,raw_lc_data, -626400, -15.632)
lc_trhust_filt = filtering(lc_lox_tank, lc_filter_order, lc_filter_cutoff)


start = 500                     # start of useful mass data intervall
end = 800                       # end of useful mass data intervall

ethanol_mass_flow, a_e, b_e, c_e = get_mass_flow(lc_ethanol_tank, lc_time, start, end)
lox_mass_folw, a_l, b_l, c_l = get_mass_flow(lc_lox_tank, lc_time, start, end)


# Plotting 
pressure_plot([ps_ethanol_tank_filt, ps_lox_tank_filt, ps_ethanol_manifold_filt, ps_lox_manifold_filt], 
            [ps_time, ps_time, ps_time, ps_time], ['ethanol tank', 'LOx tank', 'ethanol manifold', 'LOx manifold'])
        
pressure_plot([ps_chamber_1_filt, ps_chamber_2_filt], [ps_time, ps_time], ['chamber 1', 'chamber 2'])

loadcell_plot([lc_ethanol_tank_filt, lc_lox_tank_filt], [lc_time, lc_time], ['ethanol tank', 'LOx tank'], 'kg')

plot_spectrum(ps_time, ps_chamber_1, 1000, 'Spectrum Chamber PS 1')
plot_spectrum(ps_time, ps_chamber_2, 1000, 'Spectrum Chamber PS 2')