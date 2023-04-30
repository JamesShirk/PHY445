# for general math functions
import numpy as np
# for fitting the gamma spectra
from scipy.optimize import curve_fit
# for taking command arguments
import sys
# data is all in dataframes
import pandas as pd
# for plotting
import matplotlib.pyplot as plt
import matplotlib
# for intergrating peaks for # of counts
from scipy.integrate import quad
# custom modules
import fitPeak as fp
import calib
from eff import eff
from crossSectionFits import thompson, kleinNishina, fitXS
import Ne


def initialize_data(fname, names, materials):
    df = pd.read_csv(fname, header = [0, 1], skiprows=[-1])
    for i in names:
        for j in materials:
            df[i+j+"bcksub"] = df[i, j] - df[i, "bck"]

    return df

def counts(fname, names, materials):
    dat = initialize_data(fname, names, materials)
    dat = dat.loc[20:]

    # ----------
    # fit each spectrum to polynomial background + gaussian
    # ----------
    fitParams = np.empty((2, 2, 6))
    fitError = np.empty((2, 2, 6, 6))
    for i in range(len(names)):
        for j in range(len(materials)):
            elem = dat[names[i]+materials[j]+"bcksub"]
            x = np.array(elem.index.tolist())
            # hopefully height of photopeak > height of backscatter peak 
            h0 = np.max(elem)
            popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, p = h0/1.2, h=h0, peaks = True))#, sigma = 1/np.sqrt(np.abs(elem)+0.001))#, sigma=np.sqrt(np.abs(elem)+0.001))

            fitParams[i][j] = popt
            fitError[i][j] = pcov

    fig, ax = plt.subplots(2, 2, figsize = (15, 15))
    for i in range(len(names)):
        for j in range(len(materials)):
            ax[i][j].errorbar(
                dat[names[i]+materials[j]+"bcksub"].index.tolist(),
                dat[names[i]+materials[j]+"bcksub"], 
                yerr = np.sqrt(np.abs(dat[names[i]+materials[j]+"bcksub"])+0.001),
                marker = "v", ls = 'None', fillstyle = 'none', markersize = 3
            )
            ax[i][j].plot(x, fp.fit(x, fitParams[i][j]), "--")
            #ax[i][j].set_ylim([0, 1e4])
            #ax[i][j].set_yscale("symlog")
            ax[i][j].grid(visible=True)
            ax[i][j].margins(x=0)
            ax[i][j].set_xlim([0, 400])
            ax[i][j].set_title(names[i]+materials[j])
    fig.supylabel("Counts")
    fig.supxlabel("Channel")
    plt.tight_layout()
    plt.savefig("plots/diffMaterials.png", dpi = 400)
    plt.close()

    counts = np.empty((2, 2))
    countsE = np.empty((2, 2))
    Nsigma = 10

    for i in range(len(names)):
        for j in range(len(materials)):
            gauss = np.array(fitParams[i][j][3:], dtype = np.float64).reshape(3, len(fitParams[i][j][3:])//3)
            gaussE = np.sqrt(np.diagonal(fitError[i][j][3:, 3:])).reshape(3, len(fitError[i][j][3:])//3)

            # integrate from mean-nsigma*sigma to mean+nsigma*sigma
            # in general Nsigma=5 should be sufficient
            #N = quad(fp.gauss, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], args = [gauss[0][0], gauss[1][0], gauss[2][0]])
            N = quad(fp.gauss, gauss[1]-Nsigma*gauss[2], gauss[1]+Nsigma*gauss[2], args = (gauss[0], gauss[1], gauss[2]))
            N_up = quad(fp.gauss, gauss[1]-Nsigma*(gauss[2]+gaussE[2]), gauss[1]+Nsigma*(gauss[2]+gaussE[2]), args = (gauss[0]+gaussE[0], gauss[1], gauss[2]+gaussE[2]))
            N_down = quad(fp.gauss, gauss[1]-Nsigma*(gauss[2]-gaussE[2]), gauss[1]+Nsigma*(gauss[2]-gaussE[2]), args = (gauss[0]-gaussE[0], gauss[1], gauss[2]-gaussE[2]))
            counts[i][j] = N[0]
            countsE[i][j] = N_up[0]-N_down[0]

    R = np.divide(counts.T[0], counts.T[1])
    R_e = np.multiply(R, np.sqrt(np.divide(countsE.T[0], counts.T[0])**2 + np.divide(countsE.T[1], counts.T[1])**2))
    
    plt.errorbar([25, 45], R, yerr = R_e, marker = "v", ls = 'None', fillstyle = 'none', markersize = 3)
    plt.grid()
    plt.xlim([20, 50])
    plt.ylabel("Ratio of counts (Cu/Al)")
    plt.xlabel("Angle (Deg)")
    plt.savefig("plots/ratio_materials.png")


if __name__ == "__main__":
    names = ["25", "45"]
    materials = ["cu", "al"]
    counts(sys.argv[1], names, materials)