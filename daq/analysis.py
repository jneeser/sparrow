import numpy as np
import thermo
from scipy import signal, fft
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d


def filtering(data, order, cutoff):

    # gaussian blur for outlier filtering 
    data = gaussian_filter1d(data, 4)
    b, a = signal.butter(order, cutoff)
    
    # lowpass butterworth filter applied three times
    zi = signal.lfilter_zi(b, a)
    z1,_ = signal.lfilter(b, a, data, zi=zi*data[0])
    z2,_ = signal.lfilter(b, a, z1, zi=zi*z1[0])
    z3,_ = signal.lfilter(b, a, z2, zi=zi*z2[0])

    return z3


def get_ps_data(col, dataset, slope, offset):
    return dataset[:,col]*slope + offset 


def get_lc_data(col, dataset, slope, offset):
    return dataset[:,col]*slope + offset 


def get_time(col, dataset, ticks, start):
    slope = 1/ticks
    offset = dataset[:,col][0] - start*1e3
    return (dataset[:,col] - offset)*slope 


def get_dyn_pressure(mass_flow, diameter, rho):
    A = np.pi*diameter**2 / 4
    v = mass_flow / (A * rho)
    dyn_pressure = 1/2 * rho * v**2
    return dyn_pressure


def get_mass_flow(mass, time, start, end):
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    def ddx_func(x, a, b):
        return -a*b *np.exp(-b * x)

    popt,_ = curve_fit(func, time[start:end], mass[start:end])
    mass_flow_arr = - ddx_func(time[start:end], popt[0], popt[1])

    return mass_flow_arr, popt[0], popt[1], popt[2]


def get_discharge_coeff(mass_flow, area, pressure, rho):
    pressure_drop = pressure * 1e5 - 101325
    cd_arr = mass_flow / (area * np.sqrt(2 * rho * pressure_drop))
    
    return cd_arr

