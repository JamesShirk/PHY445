import glob
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(t, M0, t2p):
    return np.abs(M0*(1- 2*np.exp(-t/t2p)))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def initialize_data(dir):
    dir += "/t1*"
    files = glob.glob(dir)
    #files = [i.rsplit("/", 1)[1] for i in files]
    nums = np.array([int(i.rsplit(".", 1)[0].rsplit("_", 1)[1]) for i in files])
    idx = np.argsort(nums)
    nums = nums[idx]
    files = np.array(files)[idx]
    dfs = [pd.read_csv(i, skiprows = [-1]) for i in files]
    max = []
    for df in dfs:
        # get first row as parameters
        initParams = df.iloc[0, :].values.tolist()
        # get only 1st and second column
        df = df.iloc[1:, 0:2]
        # correct time scale
        df = df.astype({"X": float, "CH1":float})
        #df["X"] = df["X"] * float(initParams[4])
        #df["X"] = df["X"] + float(initParams[3])
        max.append(df.loc[df["CH1"].idxmax()]["CH1"])
    max = np.array(max)
    #idx = (nums-nums[np.argmin(max)] >= 0).nonzero()[0]
    #nums = nums[idx]-20
    #max = max[idx]
    x = np.linspace(nums[0], nums[-1], 1000)

    popt, pcov = curve_fit(f, nums, max, p0 = [max.max(), nums[find_nearest(max, max.min()*np.e)]])


    plt.plot(nums, max, "o")
    plt.plot(x, f(x, popt[0], popt[1]))
    plt.savefig("testt1.png")
    print(popt)

if __name__ == "__main__":
    initialize_data(sys.argv[1])