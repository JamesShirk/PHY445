import numpy as np
from scipy.optimize import curve_fit
import sys
import pandas as pd
import matplotlib.pyplot as plt

import fitPeak as fp

# def fitWrapper(x, N, *args):
#     bck = args[0]
#     gaus = []
#     for i in range(N):
#         gaus.append(args[3*i, 3*i+2])
#     return fit(x, bck, gaus)



def initialize_data(fname):
    df = pd.read_csv(fname, skiprows=[1, -1])
    # removes rows with all zeros, but preserves row numbers
    no_zero = df.loc[(df>20).any(axis=1)]
    return no_zero

def calibration(fname):
    dat = initialize_data(fname)

    calib = ['\'Fe\'', '\'Zn\'', '\'Cu\'']
    color = ["red", "green", "blue"]
    colorFit = ["purple", "orange", "m"]
    fitParams = []
    fitError = []

    for i in calib:
        elem = dat[i]
        x = np.array(elem.index.tolist())
        popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem))
        fitParams.append(popt)
        fitError.append(pcov)

    fig, ax = plt.subplots(3, 1)
    for i in range(len(calib)):
        ax[i].plot(dat[calib[i]].index.tolist(),dat[calib[i]], color = color[i])
        ax[i].plot(x, fp.fit(x, fitParams[i]), "--", color = colorFit[i])
        ax[i].set_ylim([0, 1e4])
        ax[i].set_yscale("symlog")
        ax[i].grid(visible=True)
        ax[i].margins(x=0)
        ax[i].set_xlim([200, 400])
        ax[i].set_title(calib[i])
    fig.supylabel("Counts")
    fig.supxlabel("Channel")
    plt.tight_layout()
    plt.savefig("testCalib.png", dpi = 100)
    plt.close()
    
    # in kev, order of Fe, Zn, Cu
    truePeaks = np.array([[6.40, 7.06], [8.64, 9.57], [8.05, 8.90]])
    measuredPeaks = []
    errorPeaks = []
    for i in range(len(fitParams)):
        gauss = np.array(fitParams[i][3:], dtype = np.float64).reshape(3, len(fitParams[i][3:])//3)
        #gaussE = np.array(fitError[i][3:], dtype = np.float64).reshape(3, len(fitError[i][3:])//3)
        measuredPeaks.append(gauss[1])
        errorPeaks.append(gauss[2])
    measuredPeaks = np.concatenate((measuredPeaks[0], measuredPeaks[1], measuredPeaks[2]), axis = None)
    errorPeaks = np.concatenate((errorPeaks[0], errorPeaks[1], errorPeaks[2]), axis = None)
    truePeaks = np.concatenate((truePeaks[0], truePeaks[1], truePeaks[2]), axis = None)

    p = np.argsort(measuredPeaks)

    a,b,cov_00,cov_11,cov_01,chi2 = fp.wlinear_fit(truePeaks[p],measuredPeaks[p],1.0/errorPeaks[p]**2)
    x = np.linspace(5, 11, 1000)
    y = b*x + a

    plt.errorbar(truePeaks, measuredPeaks, yerr = errorPeaks, marker = "o", ls = 'None', label = 'Measured Peaks')
    plt.plot(x, y, "--", label = 'Error Weighted Fit a:{0}, b:{1}'.format(round(a, 2), round(b, 2)))
    plt.legend()
    plt.xlabel("True Peaks (keV)")
    plt.ylabel("Channel Number")
    plt.title("Calibration with Fe, Zn, and Cu")
    plt.grid()
    plt.margins(x=0)
    plt.savefig("testCorr.png")
    


if __name__ == "__main__":
    #try:
    #initialize_data(sys.argv[1])
    calibration(sys.argv[1])
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")