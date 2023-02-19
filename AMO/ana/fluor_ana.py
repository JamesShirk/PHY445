import numpy as np
from scipy.optimize import curve_fit
import sys
import pandas as pd
import matplotlib.pyplot as plt
import fitPeak as fp
import fluor_calib as calib

def initialize_data(fname):
    df = pd.read_csv(fname, skiprows=[1, -1])
    df.columns = [col[1:-1] for col in df.columns]

    # removes rows with all zeros, but preserves row numbers
    no_zero = df.loc[(df>20).any(axis=1)]
    return df


def analysis(fname):
    dat = initialize_data(fname)
    dat = dat.iloc[25:660]

    fitParams = []
    fitError = []
    name = []

    #cols = ["Fe", "Brass", "In", "Nb"]
    #peaks = [[240, 260], [280, 305, 315, 330], [120, 560, 650], [570, 580, 600, 630]]

    for column in dat:
        elem = dat[column]
        x = np.array(elem.index.tolist())
        try:
            if column == "U" or column == "Ag":
                popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, h=10))#, peaks = peaks[i]))
            else:
                popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, h=100))#, peaks = peaks[i]))
            fitParams.append(popt)
            fitError.append(pcov)
            name.append(column)
        except:
            print(column)
        
    fig, ax = plt.subplots(len(name), 1, figsize=(20, 40))
    for i in range(len(name)):
        ax[i].plot(dat[name[i]].index.tolist(),dat[name[i]])#, color = color[i])
        ax[i].plot(x, fp.fit(x, fitParams[i]), "--")#, color = colorFit[i])
        #ax[i].set_ylim([0, 1e4])
        ax[i].set_yscale("symlog")
        ax[i].grid(visible=True)
        ax[i].margins(x=0)
        #ax[i].set_xlim([200, 400])
        ax[i].set_title(name[i])
    fig.supylabel("Counts")
    fig.supxlabel("Channel")
    plt.tight_layout()
    plt.savefig("testResult.png", dpi = 100)
    plt.close()

    a, b = calib.calibration(fname)

    measuredPeaks = []
    errorPeaks = []
    for i in range(len(fitParams)):
        gauss = np.array(fitParams[i][3:], dtype = np.float64).reshape(3, len(fitParams[i][3:])//3)
        #gaussE = np.array(fitError[i][3:], dtype = np.float64).reshape(3, len(fitError[i][3:])//3)
        measuredPeaks.append(gauss[1])
        errorPeaks.append(gauss[2])
    measuredPeaks = (np.array(measuredPeaks, dtype = object)-a)/b
    errorPeaks = (np.array(errorPeaks, dtype = object)-a)/b

    for i in range(len(measuredPeaks)):
        print(name[i])
        print(measuredPeaks[i])
        print("\n")

    #dat.plot(figsize=(20, 6))
    #plt.legend()
    #plt.savefig("test.png")

if __name__ == "__main__":
    #try:
    #initialize_data(sys.argv[1])
    analysis(sys.argv[1])
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")