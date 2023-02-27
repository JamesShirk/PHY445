import numpy as np
from scipy.optimize import curve_fit
import sys
import pandas as pd
import matplotlib.pyplot as plt
import fitPeak as fp
import fluor_calib as calib
import csv

def initialize_data(fname):
    df = pd.read_csv(fname, skiprows=[1, -1])
    df.columns = [col[1:-1] for col in df.columns]

    # removes rows with all zeros, but preserves row numbers
    no_zero = df.loc[(df>20).any(axis=1)]
    return df


def analysis(fname, plot = False):
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
    
    if plot is True:
        fig, ax = plt.subplots(len(name), 1, figsize=(20, 40))
        for i in range(len(name)):
            ax[i].plot(dat[name[i]].index.tolist(),dat[name[i]])#, color = color[i])
            ax[i].plot(x, fp.fit(x, fitParams[i]), "--")#, color = colorFit[i])
            #ax[i].set_ylim([0, 1e4])
            #ax[i].set_yscale("symlog")
            ax[i].grid(visible=True)
            ax[i].margins(x=0)
            #ax[i].set_xlim([200, 400])
            ax[i].set_title(name[i])
        fig.supylabel("Counts")
        fig.supxlabel("Channel")
        plt.tight_layout()
        plt.savefig("plots/finalRaw.png", dpi = 100)
        plt.close()

    a, b = calib.calibration(fname)

    measuredPeaks = []
    errorPeaks = []
    ampPeaks = []
    for i in range(len(fitParams)):
        gauss = np.array(fitParams[i][3:], dtype = np.float64).reshape(3, len(fitParams[i][3:])//3)
        #gaussE = np.array(fitError[i][3:], dtype = np.float64).reshape(3, len(fitError[i][3:])//3)
        ampPeaks.append(gauss[0])
        measuredPeaks.append(gauss[1])
        errorPeaks.append(gauss[2])
    measuredPeaks = (np.array(measuredPeaks, dtype = object)-a)/b
    errorPeaks = (np.array(errorPeaks, dtype = object)-a)/b

    for i in range(len(measuredPeaks)):
        print(name[i])
        print(measuredPeaks[i])
        print(errorPeaks[i])

    name_K = ["Fe", "Nb", "V", "Zn", "Zr", "Ti", "Ni", "Cu"]
    peaks_K = [[0, 1], [3, 4], [0, 1], [0, 1], [0, 1], [0, 1], [1, 2], [0, 1]]
    finalFits_K = np.empty((len(name_K), 2, 3))

    for i in range(len(name_K)):
        ind = name.index(name_K[i])
        k = 0
        for j in peaks_K[i]:
            finalFits_K[i][k] = np.array([ampPeaks[ind][j], measuredPeaks[ind][j], errorPeaks[ind][j]])
            k += 1
    name_L = ["Pb", "Au", "W", "U", "Yb", "Ag"]    
    finalFits_L = np.empty((len(name_L), 3, 3))
    peaks_L = [[1, 3, 4], [0, 5, 6], [1, 2, 5], [0, 0, 2], [2, 5, 7], [0, 0, 0]]
    for i in range(len(name_L)):
        ind = name.index(name_L[i])
        k = 0
        for j in peaks_L[i]:
            finalFits_L[i][k] = np.array([ampPeaks[ind][j], measuredPeaks[ind][j], errorPeaks[ind][j]])
            k += 1

    x = (x-a)/b
    fig, ax = plt.subplots(2, 1, figsize = (10, 10))
    for j in range(2):
        for i in range(len(finalFits_K)):
            #print(i[j][0])
            # f += i[0]*np.exp((-(x - i[1])**2) / (2*i[2]))
            ax[j].plot(x, finalFits_K[i][j][0]*np.exp((-(x - finalFits_K[i][j][1])**2) / (2*finalFits_K[i][j][2])), "-", label = name_K[i])
        ax[j].grid()
        ax[j].legend()
        ax[j].margins(x=0)
        ax[j].set_ylim(bottom=0)
    ax[0].set_title(r"$K_{\alpha}$ peaks")
    ax[1].set_title(r"$K_{\beta}$ peaks")
    fig.supxlabel("E (keV)")
    fig.supylabel("Counts")
    fig.tight_layout()
    plt.savefig("plots/final_k.png")

    fig, ax = plt.subplots(3, 1, figsize = (12, 10))
    for j in range(3):
        for i in range(len(finalFits_L)):
            if name_L[i] == "U" and j ==2:
                continue
            if name_L[i] == "Ag" and j >= 1:
                continue
            #print(i[j][0])
            # f += i[0]*np.exp((-(x - i[1])**2) / (2*i[2]))
            ax[j].plot(x, finalFits_L[i][j][0]*np.exp((-(x - finalFits_L[i][j][1])**2) / (2*finalFits_L[i][j][2])), "-", label = name_L[i])
        ax[j].grid()
        ax[j].legend()
        ax[j].margins(x=0)
        ax[j].set_ylim(bottom=0)
    ax[0].set_title(r"$L_{\alpha}$ peaks")
    ax[1].set_title(r"$L_{\beta}$ peaks")
    ax[2].set_title(r"$L_{\gamma}$ peaks")
    fig.supxlabel("E (keV)")
    fig.supylabel("Counts")
    fig.tight_layout()
    plt.savefig("plots/final_l.png")

if __name__ == "__main__":
    #try:
    #initialize_data(sys.argv[1])
    analysis(sys.argv[1], True)
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")