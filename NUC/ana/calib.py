import numpy as np
from scipy.optimize import curve_fit
import sys
import pandas as pd
import matplotlib.pyplot as plt
import fitPeak as fp
import matplotlib

def initialize_data(fname):
    df = pd.read_csv(fname, skiprows=[1, -1])
    print(df)
    print(df.columns)
    # removes quotes from name
    df.columns = [col[1:-1] for col in df.columns]
    print(df)
    return df

def calibration(fname):
    dat = initialize_data(fname)
    dat = dat.loc[100:1020]

    calib = ['Co60', 'Cs137', 'Na22', 'Ba133']

    color = ["red", "green", "blue", "orange"]
    colorFit = ["purple", "grey", "m", "brown"]
    fitParams = []
    fitError = []

    # h, prominence to find peaks good
    # in same order as calib
    init_guess = [[500, 500], [500, 500], [250, 250], [400, 150]]

    for i in range(len(calib)):
        elem = dat[calib[i]]
        x = np.array(elem.index.tolist())
        popt, pcov = curve_fit(fp.fit, xdata=x, ydata=elem, p0 = fp.initParam(elem, p = init_guess[i][1], h=init_guess[i][1]))
        fitParams.append(popt)
        fitError.append(pcov)


    if __name__ == "__main__":
        fig, ax = plt.subplots(len(calib), 1, figsize = (10, 10))
        for i in range(len(calib)):
            ax[i].plot(dat[calib[i]].index.tolist(),dat[calib[i]], color = color[i])
            ax[i].plot(x, fp.fit(x, fitParams[i]), "--", color = colorFit[i])
            #ax[i].set_ylim([0, 1e4])
            #ax[i].set_yscale("symlog")
            ax[i].grid(visible=True)
            ax[i].margins(x=0)
            #ax[i].set_xlim([200, 400])
            ax[i].set_title(calib[i])
        fig.supylabel("Counts")
        fig.supxlabel("Channel")
        plt.tight_layout()
        plt.savefig("plots/calib_fits.png", dpi = 400)
        plt.close()
    
    # in kev, order of Co60, Cs137, Na22, Ba133
    truePeaks = np.array([[1173.2, 1332.5], [661.7], [511.0, 1274.5], [302.9, 356.0]], dtype = 'object')
    
    measuredPeaks = []
    errorPeaks = []
    for i in range(len(fitParams)):
        gauss = np.array(fitParams[i][3:], dtype = np.float64).reshape(3, len(fitParams[i][3:])//3)
        #gaussE = np.array(fitError[i][3:], dtype = np.float64).reshape(3, len(fitError[i][3:])//3)
        measuredPeaks.append(gauss[1])
        errorPeaks.append(gauss[2])
    
    measuredPeaks = np.concatenate((measuredPeaks[0], measuredPeaks[1], measuredPeaks[2], measuredPeaks[3]), axis = None)
    errorPeaks = np.concatenate((errorPeaks[0], errorPeaks[1], errorPeaks[2], errorPeaks[3]), axis = None)
    truePeaks = np.concatenate((truePeaks[0], truePeaks[1], truePeaks[2], truePeaks[3]), axis = None)
    p = np.argsort(measuredPeaks)    

    a,b,cov_00,cov_11,cov_01,chi2 = fp.wlinear_fit(truePeaks[p],measuredPeaks[p],1.0/errorPeaks[p]**2)
    

    if __name__ == "__main__":
        x = np.linspace(0, 1400, 1000)
        y = b*x + a
        #print(cov_00)
        #print(cov_11)
        plt.errorbar(truePeaks, measuredPeaks, yerr = errorPeaks, marker = "o", ls = 'None', label = 'Measured Peaks')
        plt.plot(x, y, "--", 
                 label = r'Error Weighted Fit a:{0}$\pm${1}, b:{2}$\pm${3}'.format(
                    round(a, 2), 
                    round(np.sqrt(cov_00), 2), 
                    round(b, 2), 
                    round(np.sqrt(cov_11), 2)
                    )
                )
        plt.text(825, 125, r"$\chi^2$ / NDoF = {0}/{1}".format(
            round(chi2, 2),
            # Number of degrees of freedom = number of data points - number of fit parameters
            len(measuredPeaks) - 2
            ))
        plt.legend()
        plt.xlabel("True Peaks (keV)")
        plt.ylabel("Channel Number")
        plt.title("Calibration with Co60, Cs137, Ba133, and Na22")
        plt.grid()
        plt.margins(x=0)
        plt.savefig("plots/corrFinal.png")
    return a, b
    


if __name__ == "__main__":
    #try:
    #matplotlib.rcParams.update({'font.size': 18})
    #initialize_data(sys.argv[1])
    calibration(sys.argv[1])
    #except:
    #    print("Run as \'python3 fluor_calib.py datafile\'")