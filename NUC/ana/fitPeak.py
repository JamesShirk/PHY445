import numpy as np
from scipy.signal import find_peaks
from itertools import chain
import pandas as pd


def initParam(elem, p = 30, h = 120, peaks = None):
    if peaks is None:
        peaks, _ = find_peaks(elem, height = h, prominence = p)
        amps = elem.iloc[peaks].tolist()
    else:
        peaks, _ = find_peaks(elem, height = h, prominence = p)
        amps = elem.loc[peaks].tolist()
        print(amps)
    means = elem.index[peaks]
    sigmas = np.ones_like(amps)*5
    bck = [10, 0, 0]
    return np.array(list(chain(bck, amps, means, sigmas)), dtype = np.float64)

def fit(x, *args):
    # alright fine ill finally break down how this function works since i cant even remember
    # scipy curve_fit only works with a fixed number of fit parameters, but our spectra has multiple peaks
    # we fit these peaks with gaussians 3 parameters each, so we have to do some *args fuckery to get it to work

    # anyways, the paramters end up being in an n*3 + 3 array where n = number of peaks
    # array goes like [polynomial coefficients (a + bx + cx^2), amplitude_i, mean_i, stdev_i, amplitude_i+1, ...]


    if type(args[0]) is np.ndarray:
        args = np.array(args[0], dtype = np.float64)


    f = np.zeros_like(x, dtype = np.float64)
    gauss = np.array(args[3:], dtype = np.float64).reshape(3, len(args[3:])//3).T
    for i in gauss:
        f += i[0]*np.exp((-(x - i[1])**2) / (2*i[2]**2))
    #return f + args[0]*np.exp(-args[1]*x +args[2])
    return f + args[0] + args[1]*x + args[2]*x**2

# taken from stackexchange
def wlinear_fit (x,y,w) :
    """
    Fit (x,y,w) to a linear function, using exact formulae for weighted linear
    regression. This code was translated from the GNU Scientific Library (GSL),
    it is an exact copy of the function gsl_fit_wlinear.
    """
    # compute the weighted means and weighted deviations from the means
    # wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i)
    W = np.sum(w)
    wm_x = np.average(x,weights=w)
    wm_y = np.average(y,weights=w)
    dx = x-wm_x
    dy = y-wm_y
    wm_dx2 = np.average(dx**2,weights=w)
    wm_dxdy = np.average(dx*dy,weights=w)
    # In terms of y = a + b x
    b = wm_dxdy / wm_dx2
    a = wm_y - wm_x*b
    cov_00 = (1.0/W) * (1.0 + wm_x**2/wm_dx2)
    cov_11 = 1.0 / (W*wm_dx2)
    cov_01 = -wm_x / (W*wm_dx2)
    # Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2
    chi2 = np.sum (w * (y-(a+b*x))**2)
    return a,b,cov_00,cov_11,cov_01,chi2