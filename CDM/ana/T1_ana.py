import glob
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(t, M0, t1, a):
    return np.abs(M0*(1- 2*np.exp(-t/t1))) + a

# finds nearest value in array to given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def initialize_data(dir, outname):
    dir += "/t1*"
    # read in files
    files = glob.glob(dir)
    '''
    filenames have convention of t1_substance_delay.csv

    following lines get delay and sort files by delay
    '''
    nums = np.array([int(i.rsplit(".", 1)[0].rsplit("_", 1)[1]) for i in files])
    idx = np.argsort(nums)
    nums = nums[idx]
    files = np.array(files)[idx]
    dfs = [pd.read_csv(i, skiprows = [-1]) for i in files]
    max = []
    # for each file, find the peak on channel 1(i.e. amplitude after pi then pi/2 pulse)
    for df in dfs:
        # get first row as parameters
        initParams = df.iloc[0, :].values.tolist()
        # get only 1st and second column
        df = df.iloc[1:, 0:2]
        # correct time scale
        df = df.astype({"X": float, "CH1":float})
        # these lines correct time scale, doesnt actually matter for this measurement
        #df["X"] = df["X"] * float(initParams[4])
        #df["X"] = df["X"] + float(initParams[3])
        max.append(df.loc[df["CH1"].idxmax()]["CH1"])
    max = np.array(max)
    x = np.linspace(nums[0], nums[-1], 1000)
    T1_guess = nums[find_nearest(max, max.min()*np.e)]
    if outname == "wt":
        a = 1.5
        T1_guess = 3000
    else:
        a=0
    # fit w/ initial amplitude guess as max of array, initial T1 = time for signal to increase by factor of e
    popt, pcov = curve_fit(f, nums, max, p0 = [max.max(), T1_guess, a])

    # plotting
    fig, ax = plt.subplots()
    plt.plot(nums, max, "o", label = "data")
    plt.plot(x, f(x, *popt), label = "fit")
    plt.grid()
    plt.legend()
    plt.xlabel(r"$\pi-\pi/2$ signal delay (ms)")
    plt.ylabel("Amplified Signal (V)")
    plt.text(0.6, 0.7, r"$T_1 = {0} \pm {1}\quad (ms)$".format(
        #round(popt[1], 3),
        round(popt[1]//10 * 10),
        #round(np.sqrt(pcov[1][1]), 1)
        round(np.sqrt(pcov[1][1])//10 * 10)
        ), transform = ax.transAxes)
    plt.savefig("T1/main_T1_"+outname+".png")

if __name__ == "__main__":
    initialize_data(sys.argv[1], sys.argv[2])