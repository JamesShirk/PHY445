import numpy as np
from scipy.optimize import curve_fit
import sys
import pandas as pd
import matplotlib.pyplot as plt
import fitPeak as fp
import matplotlib


def initialize_data(fname, names):
    df = pd.read_csv(fname, header = [0, 1], skiprows=[-1])
    #df.columns = [col[1:-1] for col in df.columns]
    #df = df.loc[0:1020]
    for i in names:
        df[i+"bcksub"] = df[i, "main"] - df[i, "bck"]
    df1 = df.loc[:, names[0]+"bcksub":names[-1]+"bcksub"]
    return df1

def calibration(fname):
    names = ["20", "40", "90"]
    dat = initialize_data(fname, names)

    fig, ax = plt.subplots(len(names), 1, figsize = (10, 10))
    for i in range(len(names)):
        ax[i].plot(dat[names[i]+"bcksub"].index.tolist(),dat[names[i]+"bcksub"])#, color = color[i])
        #ax[i].set_ylim([0, 1e4])
        #ax[i].set_yscale("symlog")
        ax[i].grid(visible=True)
        ax[i].margins(x=0)
        #ax[i].set_xlim([200, 400])
        ax[i].set_title(names[i])
    fig.supylabel("Counts")
    fig.supxlabel("Channel")
    plt.tight_layout()
    plt.savefig("plots/main_raw.png", dpi = 400)
    plt.close()

if __name__ == "__main__":
    #try:
    #matplotlib.rcParams.update({'font.size': 18})
    #initialize_data(sys.argv[1])
    calibration(sys.argv[1])
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")