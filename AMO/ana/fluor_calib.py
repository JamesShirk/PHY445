import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.fft import fft, ifft, fftfreq
import sys
import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain

# def fitWrapper(x, N, *args):
#     bck = args[0]
#     gaus = []
#     for i in range(N):
#         gaus.append(args[3*i, 3*i+2])
#     return fit(x, bck, gaus)

def initParam(elem, peaks):
    amps = elem.iloc[peaks].tolist()
    means = elem.index[peaks]
    print(means)
    sigmas = np.ones_like(amps)*5
    bck = [10, 0, 0]
    return np.array(list(chain(bck, amps, means, sigmas)), dtype = np.float64)

def fit(x, *args):
    if type(args[0]) is np.ndarray:
        args = np.array(args[0], dtype = np.float64)
    f = np.zeros_like(x, dtype = np.float64)
    gauss = np.array(args[3:], dtype = np.float64).reshape(3, len(args[1:])//3).T
    for i in gauss:
        f += i[0]*np.exp((-(x - i[1])**2) / (2*i[2]))
    return f + args[0] + args[1]*x + args[2]*x**2

def initialize_data(fname):
    df = pd.read_csv(fname, skiprows=[1, -1])
    # removes rows with all zeros, but preserves row numbers
    no_zero = df.loc[(df>20).any(axis=1)]
    fig, ax = plt.subplots()
    #ax.plot(no_zero)#, x = "channel", y = "counts")
    ax = no_zero.plot(figsize= (20, 6))
    #fig.set_size_inches(10, 6, forward = True)
    ax.set_xlim([50, 900])
    ax.set_ylim([10, 1e4])
    ax.set_yscale("symlog")
    ax.legend()
    plt.savefig("test.png", dpi =100)
    plt.close()

    dx = 1
    
    
    peak_subset = no_zero['\'Rh\'']#.loc[200:300]


    #curve_fit(fit, no_zero.index, no_zero['\'Fe\''])
    peaks, _ = find_peaks(peak_subset, height = 10, prominence = dx*50)
    x = np.array(peak_subset.index.tolist())
    
    popt, pcov = curve_fit(fit, xdata=x, ydata=peak_subset, p0 = initParam(peak_subset, peaks))


    fig, ax = plt.subplots()
    ax = peak_subset.plot(figsize= (20, 6))
    ax.plot(x, fit(x, popt), "--")
    #ax.set_xlim([50, 300])
    ax.set_ylim([0, 1e4])
    ax.set_yscale("symlog")
    plt.savefig("testRh.png", dpi = 100)
    plt.close()

    return

    Ft = fft(no_zero['\'Fe\''])
    N = no_zero['\'Fe\''].shape[0]
    xf = fftfreq(N, 1)#[:N//2]
    print(xf)
    xf = np.where(np.abs(xf) > 0.2)
    print(xf)
    cut = len(xf)//2
    Ift = ifft(Ft[xf])
    #xf = np.argmax(np.abs(xf)>0.2)

    fig, ax = plt.subplots()
    ax.plot(xf[0], Ift)
    #ax.plot(xf, 2.0/N * np.abs(Ft[0:N//2].imag), "--")
    plt.savefig("testFourier.png")



if __name__ == "__main__":
    #try:
    initialize_data(sys.argv[1])
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")