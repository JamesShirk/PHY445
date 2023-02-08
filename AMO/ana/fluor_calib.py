import numpy as np
import scipy as sci
import sys
import pandas as pd
import matplotlib.pyplot as plt

def initialize_data(fname):
    df = pd.read_csv(fname, skiprows=[1, -1])
    # removes rows with all zeros, but preserves row numbers
    no_zero = df.loc[(df!=0).any(axis=1)]
    fig, ax = plt.subplots()
    #ax.plot(no_zero)#, x = "channel", y = "counts")
    ax = no_zero.plot(figsize= (20, 6))
    #fig.set_size_inches(10, 6, forward = True)
    ax.set_xlim([50, 900])
    ax.set_ylim([10, 1e4])
    ax.set_yscale("symlog")
    ax.legend()
    plt.savefig("test.png", dpi =100)
    


if __name__ == "__main__":
    #try:
    initialize_data(sys.argv[1])
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")