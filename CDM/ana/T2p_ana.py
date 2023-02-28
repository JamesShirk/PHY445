import pandas as pd
import matplotlib.pyplot as plt
import sys
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import numpy as np


def f(t, A, T2p, B):
    return A*np.exp(-t/T2p)+B

def initialize_data(fname):
    # read in csv
    df = pd.read_csv(fname, skiprows = [-1])
    # get first row as parameters
    initParams = df.iloc[0, :].values.tolist()
    # get only 1st and second column
    df = df.iloc[1:, 0:2]
    # correct time scale
    df = df.astype({"X": float, "CH1":float})
    df["X"] = df["X"] * float(initParams[3])
    df["X"] = df["X"] + float(initParams[2])
    return df

def t2prime(fname, outname):
    dat = initialize_data(fname)
    peaks, _ = find_peaks(dat["CH1"], height = 1.0)
    
    plt.plot(dat["X"], dat["CH1"])
    plt.grid()
    plt.xlabel("Time (ms)")
    plt.ylabel("Amplified Signal (V)")
    plt.savefig("./T2Prime/raw_t2p_"+outname+".png")
    plt.close()
    
    onlyPeaks = dat.iloc[peaks[1:]]
    x = np.array(onlyPeaks["X"].tolist())
    popt, pcov = curve_fit(f, onlyPeaks["X"], onlyPeaks["CH1"], p0 = [onlyPeaks.iloc[0]["CH1"], x[(len(x)-1)//2], 1])
    x_fit = np.linspace(x[0], x[-1], 100)
    
    fig, ax = plt.subplots()
    plt.plot(onlyPeaks["X"], onlyPeaks["CH1"], "o", label = "Data")
    plt.plot(x_fit, f(x_fit, popt[0], popt[1], popt[2]), label = "fit")
    plt.text(0.5, 0.5, r"$T_2^\prime = {0:.2f} \pm {1:.2f}$ (ms)".format(
        round(popt[1]*1e3, 2), 
        round(np.sqrt(pcov[1][1])*1e3, 2)
        ), transform = ax.transAxes)
    plt.legend()
    plt.grid()
    plt.xlabel("Time (ms)")
    plt.ylabel("Amplified Signal (V)")
    plt.savefig("./T2Prime/main_t2p_"+outname+".png")
    plt.close()
    
    fig, ax = plt.subplots()
    plt.plot(dat["X"], dat["CH1"], label = "All data")
    plt.plot(onlyPeaks["X"], onlyPeaks["CH1"], "o", label = "Extracted peaks")
    plt.plot(x_fit, f(x_fit, popt[0], popt[1], popt[2]), label = "Fit")
    plt.grid()
    plt.legend()
    plt.text(0.6, 0.7, r"$T_2^\prime =$ {} $ms \pm$ {}".format(round(popt[1]*1e3, 2), round(np.sqrt(pcov[1][1])*1e3, 2)), transform = ax.transAxes)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplified Signal (V)")
    plt.savefig("./T2Prime/tot_t2_"+outname+".png")
    plt.close()


if __name__ == "__main__":
    t2prime(sys.argv[1], sys.argv[2])