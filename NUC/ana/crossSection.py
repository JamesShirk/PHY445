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

def bckOnly(x, a, b, c):
    return a + b*x + c*x**2

def initialize_data(fname, names):
    df = pd.read_csv(fname, header = [0, 1], skiprows=[-1])
    #df.columns = [col[1:-1] for col in df.columns]
    #df = df.loc[0:1020]
    for i in names:
        df[i+"bcksub"] = df[i, "main"] - df[i, "bck"]
    df1 = df.loc[:, names[0]+"bcksub":names[-1]+"bcksub"]
    return df1

def Ngamma(runTime = 900):
    df2 = pd.read_csv("../data/0deg.out", header = [0], skiprows=[1, -1])
    x = np.array(df2["0"].index.tolist())
    h0 = 10000
    popt, pcov = curve_fit(fp.fit, xdata=x, ydata=df2["0"], p0 = fp.initParam(df2["0"], p = h0, h=h0, peaks=True))
    gauss = np.array(popt[3:], dtype = np.float64).reshape(3, len(popt[3:])//3)
    Nsigma = 10
    N = quad(fp.fit, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], popt)
    bck = quad(bckOnly, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], (popt[0], popt[1], popt[2]))
    return (N[0] - bck[0])/ runTime

def counts(fname):
    # list of angles to calculate using
    names = ["20", "30", "40", "50", "90"]
    runTime = np.array([900, 900, 900, 1200, 1200])
    dat = initialize_data(fname, names)
    

    # ----------
    # fit each spectrum to polynomial background + gaussian
    # ----------
    fitParams = []
    fitError = []
    for i in range(len(names)):
        elem = dat[names[i]+"bcksub"]
        x = np.array(elem.index.tolist())
        # hopefully height of photopeak > height of backscatter peak 
        h0 = np.max(elem)-1
        popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, p = h0/1.5, h=h0))#, sigma=np.sqrt(np.abs(elem)+0.001))
        fitParams.append(popt)
        fitError.append(pcov)

    # ----------
    # plot raw spectra + fit
    # ----------

    fig, ax = plt.subplots(len(names), 1, figsize = (10, 10))
    for i in range(len(names)):
        ax[i].errorbar(
            dat[names[i]+"bcksub"].index.tolist(),
            dat[names[i]+"bcksub"], 
            yerr = np.sqrt(np.abs(dat[names[i]+"bcksub"])+0.001),
            marker = "v", ls = 'None', fillstyle = 'none', markersize = 3
        )
        ax[i].plot(x, fp.fit(x, fitParams[i]), "--")
        #ax[i].set_ylim([0, 1e4])
        #ax[i].set_yscale("symlog")
        ax[i].grid(visible=True)
        ax[i].margins(x=0)
        ax[i].set_xlim([0, 400])
        ax[i].set_title(names[i])
    fig.supylabel("Counts")
    fig.supxlabel("Channel")
    plt.tight_layout()
    plt.savefig("plots/main_raw.png", dpi = 400)
    plt.close()

    # ----------
    # calculate number of counts
    # ----------

    counts = []
    measuredPeaks = []
    errorPeaks = []
    Nsigma = 10


    for i in range(len(fitParams)):
        gauss = np.array(fitParams[i][3:], dtype = np.float64).reshape(3, len(fitParams[i][3:])//3)
        # integrate from mean-nsigma*sigma to mean+nsigma*sigma
        # in general Nsigma=5 should be sufficient
        N = quad(fp.fit, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], fitParams[i])
        bck = quad(bckOnly, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], (fitParams[i][0], fitParams[i][1], fitParams[i][2]))
        counts.append(N[0]-bck[0])
        measuredPeaks.append(gauss[1])
        errorPeaks.append(gauss[2])

    # returns counts per second since each data point wasnt run for equal amount of time
    # runTime is defined at start of function

    return np.divide(np.array(counts), runTime), measuredPeaks, errorPeaks


def kleinNishina(theta,E):
    g = E/511
    p1 = (1 + np.cos(theta * np.pi/180)**2)/((1+g*(1-np.cos(theta * np.pi/180)))**2)
    p2 = (1 + (g**2 * (1-np.cos(theta * np.pi/180))**2)/ ((1+np.cos(theta * np.pi/180)**2)*(1+g*(1-np.cos(theta * np.pi/180)))))
    return (2.82e-13**2/2)*p1*p2 * 1e27

def thompson(theta):
    return (2.82e-13**2/2)*(1+np.cos(theta * np.pi/180)**2)* 1e27

def crossSection(counts, theta, ngamma, nelectron):

    dsigma = counts* 1e27 / (ngamma*nelectron*0.03)
    x = np.linspace(theta[0], theta[-1], 100)
    kn = kleinNishina(x, 662)
    th = thompson(x)
    dsigma *= kn[0]/dsigma[0]
    plt.plot(theta, dsigma, ls = 'None', marker = 'o', label = "Measured")
    plt.plot(x, kn, "--", label = "Klein-Nishina")
    plt.plot(x, th, ls = "dashdot", label = "Thompson")
    plt.yscale("log")
    plt.legend()
    plt.xlabel(r"$\theta$ (deg)")
    plt.ylabel(r"$d\sigma/d\Omega$ (mb)")
    plt.savefig("plots/differentialCrossSection.png", dpi = 300)
    plt.close()


def energy(peaks, peaksError, runs):
    a, b, aerr, berr  =  calib.calibration("../data/calib.out")
    peaks = (peaks - a)/(b*1e3)
    peaksError = (1/b) * np.sqrt(peaksError**2 + aerr + berr*peaks**2) / (1e3)
    #print(peaksError)

    #peaksError = (peaksError - a)/(b*1e3)

    #print(peaks)

    x = 1-np.cos(runs*np.pi/180)
    y = (1/peaks - 1/662).T[0]
    yerr = np.divide(peaksError, peaks).T[0]

    a1, b1, cov_00,cov_11,cov_01,chi2 = fp.wlinear_fit(x, y, 1/yerr**2)
    xFit = np.linspace(x[0], x[-1], 100)

    # propagate error given x = y-a/b where all of y, a, b have error

    plt.errorbar(x, y, yerr = yerr, ls = 'None', marker = "o", label = "Measured")
    plt.plot(xFit, a1 + b1*xFit, "--", label = "Fit")
    plt.ylabel(r"1/E (MeV$^{-1}$)")
    plt.xlabel(r"1-$\cos(\theta)$")
    plt.legend()
    plt.text(x = 0.63, y = 2.2, s = r"$m_e$ = {0} $\pm$ {1} KeV".format(
        round(1/b1 * 1e3, 0),
        round(np.sqrt(cov_11)/b1 * 1e3, 0)
    ))
    plt.text(x=0.63, y=1.8, s = r"$\chi^2$/NDF = {0}/{1}".format(
        round(chi2, 2),
        len(runs) - 2
    ))
    plt.grid()
    plt.tight_layout()
    plt.savefig("plots/electronMass.png", dpi = 300)
    plt.close()


if __name__ == "__main__":
    #try:
    #matplotlib.rcParams.update({'font.size': 18})
    runs = np.array([20, 30, 40, 50, 90])
    N, peaks, peaksError = counts(sys.argv[1])
    #Ngamma()
    #crossSection(N, runs, Ngamma(), 4e26)
    energy(np.array(peaks), np.array(peaksError), runs)