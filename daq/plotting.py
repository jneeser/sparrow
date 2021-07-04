from matplotlib import pyplot as plt
import numpy as np


def plot_spectrum(time, data, sampling_freq, title):
    NFFT = 1024
    Fs = sampling_freq


    fig, ax = plt.subplots(nrows=1)
    Pxx, freqs, bins, im = ax.specgram(data, NFFT=NFFT, Fs=Fs, noverlap=900)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('frequency [Hz]')

    plt.title(title)
    plt.show()  


def pressure_plot(data, times, labels):

    for i in range(len(data)):
        plt.plot(times[i], data[i], label = labels[i])

    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('pressure [bar]')
    plt.title('pressure data')
    plt.grid()
    plt.show()    


def loadcell_plot(data, times, labels, unit):

    for i in range(len(data)):
        plt.plot(times[i], data[i], label = labels[i])

    plt.legend(loc='best')
    plt.xlabel('time [s]')

    if unit == 'kg':
        plt.ylabel('mass [kg]')

    elif unit == 'N':
        plt.ylabel('Force [N]')

    plt.title('loadcell data')
    plt.grid()
    plt.show()


def massflow_plot(data, times, coefficients, labels):
    
    for i in range(len(data)):
        plt.plot(data[i], abs(-coefficients[i][0]*coefficients[i][1] *np.exp(-coefficients[i][1] * times[i])), label = labels[i])

    plt.xlabel('tank pressure [bar]')
    plt.ylabel('mass flow [kg/s]')
    plt.title('injector mass flow')
    plt.grid()
    plt.show()