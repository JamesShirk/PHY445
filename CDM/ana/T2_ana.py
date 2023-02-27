import pandas as pd
import matplotlib.pyplot as plt
import sys
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import numpy as np


def f(t, tp, A, T2, B):
    return A*np.exp(-(t-tp)/T2) + B

def initialize_data(fname):
    # read in csv
    df = pd.read_csv(fname, skiprows = [-1])
    # get first row as parameters
    initParams = df.iloc[0, :].values.tolist()
    # get only 1st and second column
    df = df.iloc[1:, 0:2]
    # correct time scale
    df = df.astype({"X": float, "CH1":float})
    df["X"] = df["X"] * float(initParams[4])
    df["X"] = df["X"] + float(initParams[3])

    return df.loc[df["X"] > 0]

def get_closest(df, val, col):
    return(df[col]-val).abs().idxmin()


def t2prime(fname):
    dat = initialize_data(fname)
    # Should only find one peak
    peak = find_peaks(dat["CH1"], height=1, prominence = 3, distance = 100)[0]
    if len(peak) > 1:
        print("Error, found multiple peaks, check data")
        return
    dat = dat.iloc[peak[0]+10:peak[0]+200]
    A = dat.iloc[0]["CH1"]
    t0 = dat.iloc[0]["X"]
    #initial guess for t2 is time after which amplitude drops by e
    t2_g = float(dat.loc[dat.index == get_closest(dat, A/np.e, "CH1")]["X"])
    popt, pcov = curve_fit(f, dat["X"], dat["CH1"], p0 = [t0, A, t2_g, 0])
    dat.plot("X", "CH1", marker="o", fillstyle = 'none', markersize=2)
    print(popt[2]*1e3)
    plt.plot(dat["X"], f(dat["X"], popt[0], popt[1], popt[2], popt[3]))
    plt.savefig("testT2.png")


if __name__ == "__main__":
    t2prime(sys.argv[1])