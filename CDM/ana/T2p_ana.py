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

def t2prime(fname):
    dat = initialize_data(fname)
    peaks, _ = find_peaks(dat["CH1"], height = 1.0)
    onlyPeaks = dat.iloc[peaks[1:]]
    x = np.array(onlyPeaks["X"].tolist())
    popt, pcov = curve_fit(f, onlyPeaks["X"], onlyPeaks["CH1"], p0 = [onlyPeaks.iloc[0]["CH1"], x[(len(x)-1)//2], 1])
    plt.plot(onlyPeaks["X"], onlyPeaks["CH1"], "o")
    x_fit = np.linspace(x[0], x[-1], 100)
    plt.plot(x_fit, f(x_fit, popt[0], popt[1], popt[2]))
    plt.savefig("testt2p.png")
    print(popt)


if __name__ == "__main__":
    t2prime(sys.argv[1])