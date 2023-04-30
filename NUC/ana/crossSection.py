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


def initialize_data(fname, names):
    df = pd.read_csv(fname, header = [0, 1], skiprows=[-1])
    for i in names:
        df[i+"bcksub"] = df[i, "main"] - df[i, "bck"]

    return df

def Ngamma(runTime = 900):
    df2 = pd.read_csv("../data/0deg.out", header = [0], skiprows=[1, -1])
    x = np.array(df2["0"].index.tolist())
    h0 = 10000
    popt, pcov = curve_fit(fp.fit, xdata=x, ydata=df2["0"], p0 = fp.initParam(df2["0"], p = h0, h=h0, peaks=None))
    gauss = np.array(popt[3:], dtype = np.float64)#.reshape(3, len(popt[3:])//3)
    Nsigma = 10
    #N = quad(fp.gauss, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], args = [gauss[0][0], gauss[1][0], gauss[2][0]])
    N = quad(fp.gauss, gauss[1]-Nsigma*gauss[2], gauss[1]+Nsigma*gauss[2], args = (gauss[0], gauss[1], gauss[2]))
    #bck = quad(bckOnly, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], (popt[0], popt[1], popt[2]))
    #return (N[0] - bck[0])/ runTime
    #return N[0]/runTime
    return N[0]

def counts(fname, runs, runTime):
    # list of angles to calculate using
    names = [str(i) for i in runs]
    dat = initialize_data(fname, names)
    dat = dat.loc[20:]

    # ----------
    # fit each spectrum to polynomial background + gaussian
    # ----------
    fitParams = []
    fitError = []
    for i in range(len(names)):
        elem = dat[names[i]+"bcksub"]
        x = np.array(elem.index.tolist())
        # hopefully height of photopeak > height of backscatter peak 
        h0 = np.max(elem)
        if "120" in names[i] or "130" in names[i]:
            #print('hi')
            popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, p = h0/1.2, h=h0, peaks = True), sigma = 1/np.log(np.sqrt(np.abs(elem)+0.001)))#, sigma=np.sqrt(np.abs(elem)+0.001))
        else:
            popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, p = h0/1.2, h=h0, peaks = True))#, sigma = 1/np.sqrt(np.abs(elem)+0.001))#, sigma=np.sqrt(np.abs(elem)+0.001))

        fitParams.append(popt)
        fitError.append(pcov)

    # ----------
    # plot raw spectra + fit
    # ----------

    fig, ax = plt.subplots()
    ax.errorbar(
        dat["20bcksub"].index.tolist(),
        dat["20bcksub"], 
        yerr = np.sqrt(np.abs(dat["20bcksub"])+0.001),
        marker = "v", ls = 'None', fillstyle = 'none', markersize = 3, label = "Background subtracted"
    )
    ax.plot(dat["20", "main"].index.tolist(), dat["20", "main"], "--", label = "Scatterer")
    ax.plot(dat["20", "bck"].index.tolist(), dat["20", "bck"], label = "No Scatterer")
    ax.plot(dat["20bcksub"].index.tolist(), fp.gauss(dat["20bcksub"].index.tolist(), fitParams[0][3], fitParams[0][4], fitParams[0][5]), label = "Main peak")
    ax.grid(visible=True)
    ax.margins(x=0)
    ax.set_xlim([0, 400])
    ax.set_title("20 degree")
    ax.legend()
    ax.set_ylabel("Counts")
    ax.set_xlabel("Channel")
    plt.tight_layout()
    plt.savefig("plots/20deg_raw.png", dpi = 400)
    plt.close()



    fig, ax = plt.subplots(len(names), 1, figsize = (10, 30))
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
    countsE = []
    measuredPeaks = []
    errorPeaks = []
    Nsigma = 10




    for i in range(len(fitParams)):
        gauss = np.array(fitParams[i][3:], dtype = np.float64).reshape(3, len(fitParams[i][3:])//3)
        gaussE = np.sqrt(np.diagonal(fitError[i][3:, 3:])).reshape(3, len(fitError[i][3:])//3)

        # integrate from mean-nsigma*sigma to mean+nsigma*sigma
        # in general Nsigma=5 should be sufficient
        #N = quad(fp.gauss, gauss[1][0]-Nsigma*gauss[2][0], gauss[1][0]+Nsigma*gauss[2][0], args = [gauss[0][0], gauss[1][0], gauss[2][0]])
        N = quad(fp.gauss, gauss[1]-Nsigma*gauss[2], gauss[1]+Nsigma*gauss[2], args = (gauss[0], gauss[1], gauss[2]))
        N_up = quad(fp.gauss, gauss[1]-Nsigma*(gauss[2]+gaussE[2]), gauss[1]+Nsigma*(gauss[2]+gaussE[2]), args = (gauss[0]+gaussE[0], gauss[1], gauss[2]+gaussE[2]))
        N_down = quad(fp.gauss, gauss[1]-Nsigma*(gauss[2]-gaussE[2]), gauss[1]+Nsigma*(gauss[2]-gaussE[2]), args = (gauss[0]-gaussE[0], gauss[1], gauss[2]-gaussE[2]))
        counts.append(N[0])
        countsE.append(N_up[0]-N_down[0])
        measuredPeaks.append(gauss[1])
        errorPeaks.append(gauss[2])

    # returns counts per second since each data point wasnt run for equal amount of time
    # runTime is defined at start of function

    #print(counts)
    return np.divide(np.array(counts), runTime)*900, np.divide(np.array(countsE), runTime)*900,measuredPeaks, errorPeaks


def crossSection(counts, countsError, theta, ngamma, nelectron, peaks):

    efficiency, peak_to_total = eff()
    intercept, slope, aerr, berr  =  calib.calibration("../data/calib.out")
    peaks = (peaks - intercept)/(slope)

    dsigmaE = np.divide(countsError, counts)


    dsigma = counts / (ngamma*nelectron)


    for i in range(len(dsigma)):
       dsigma[i] = dsigma[i]/efficiency(peaks[i])
       dsigma[i] = dsigma[i]/peak_to_total(peaks[i])



    dsigma *= 1e24



    #x = np.linspace(theta[0], theta[-1], 100)
    x = np.linspace(0, 180, 100)
    kn = kleinNishina(x, 662) * 1e28
    th = thompson(x) * 1e28


    popt, pcov = curve_fit(lambda x, a, b: fitXS(x, a, E = 662, b=b), theta, dsigma, p0 = [1, 0.1])

    dsigmaE = np.multiply(dsigma, dsigmaE)

    plt.errorbar(theta, dsigma, yerr = dsigmaE, ls = 'None', marker = 'o', label = "Measured")
    plt.plot(x, kn, "--", label = "Klein-Nishina")
    plt.plot(x, th, ls = "dashdot", label = "Thompson")
    plt.plot(x, fitXS(x, popt[0], E=662, b = popt[1]), ls = "--", label = "Fit")
    plt.grid()
    plt.legend()
    plt.xlabel(r"$\theta$ (deg)")
    plt.ylabel(r"$d\sigma/d\Omega$ (b)")
    plt.tight_layout()
    plt.savefig("plots/differentialCrossSection.png", dpi = 200)
    plt.close()

    return np.sqrt(popt[0] * 2/1e28), np.sqrt(np.sqrt(pcov[0][0]) * 2/1e28)

def energy(peaks, peaksError, runs):
    # peaks = peaks
    # peaksError = peaksError
    # runs = runs
    a, b, aerr, berr  =  calib.calibration("../data/calib.out")
    peaks = (peaks - a)/(b*1e3)
    peaksError = (1/b) * np.sqrt(peaksError**2 + aerr + berr*peaks**2) / (1e3)

    x = (1-np.cos(runs*np.pi/180))
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

    return 1/b1 * 1e3, np.sqrt(cov_11)/b1 * 1e3

if __name__ == "__main__":
    #try:
    #matplotlib.rcParams.update({'font.size': 18})
    runTime = np.array([900, 1200, 900, 1200, 900, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200])
    runs = np.array([    20,   25,  30,   35,  40,   50,   60,   70,   80,   90,  100, 110,  120,  130])
    N, Nerror, peaks, peaksError = counts(sys.argv[1], runs, runTime)
    #Ngamma()
    # in cm i think lol maybe
    area = (1.5**2*np.pi) * 6.45
    r_o, r_o_err = crossSection(N, Nerror, runs, Ngamma()/area, Ne.electrons(), peaks)
    # in keV
    m_e, m_e_err = energy(np.array(peaks), np.array(peaksError), runs)
    m_e *= 1e3
    m_e_err *= 1e3

    print(m_e)
    print(m_e_err)
    print(r_o)
    print(r_o_err)

    # in eV*s
    hbar = 6.582e-16
    # in m/s
    c = 2.99e8
    print((m_e*r_o/(hbar*c)))
    print((1/(hbar*c))*np.sqrt(r_o**2 * m_e_err**2 + m_e**2 * r_o_err**2))
    print("Fine structure constant is: 1/{} +- 1/{}".format(
        round(1/(m_e*r_o/(hbar*c)), 1),
        round(1/ ((1/(hbar*c))*np.sqrt(r_o**2 * m_e_err**2 + m_e**2 * r_o_err**2)), 1)
    ))
    print(1/ ((m_e*r_o/(hbar*c)) + (1/(hbar*c))*np.sqrt(r_o**2 * m_e_err**2 + m_e**2 * r_o_err**2)))
    print(1/ ((m_e*r_o/(hbar*c)) - (1/(hbar*c))*np.sqrt(r_o**2 * m_e_err**2 + m_e**2 * r_o_err**2)))
