import numpy as np
from scipy.signal import find_peaks
from itertools import chain
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


def initParam(elem, dx = 1, h = 120, peaks = None):
    if peaks is None:
        peaks, _ = find_peaks(elem, height = h, prominence = dx*30)
        amps = elem.iloc[peaks].tolist()
    else:
        amps = elem.loc[peaks].tolist()
        print(amps)
    means = elem.index[peaks]
    sigmas = np.ones_like(amps)*5
    bck = [10, 0, 0]
    return np.array(list(chain(bck, amps, means, sigmas)), dtype = np.float64)

def fit(x, *args):
    if type(args[0]) is np.ndarray:
        args = np.array(args[0], dtype = np.float64)
    f = np.zeros_like(x, dtype = np.float64)
    gauss = np.array(args[3:], dtype = np.float64).reshape(3, len(args[3:])//3).T
    for i in gauss:
        f += i[0]*np.exp((-(x - i[1])**2) / (2*i[2]))
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



x = np.array([0.008514,0.009933,0.011919,0.014899]);
y = np.array([0.057,0.069,0.081,0.101]);
delta = np.array([0.000003,0.000004,0.000005,0.000008]);
w = np.array([1/(0.000003)**2, 1/(0.000004)**2, 1/(0.000005)**2, 1/(0.000008)**2 ])
deltax = np.array([0.007,0.007,0.004,0.005]);
wx = np.array([1/(0.007)**2, 1/(0.007)**2, 1/(0.004)**2, 1/(0.005)**2 ])

# 

#plt.plot(x,y);
#print(wlinear_fit(x,y,w))
#intercept,slope,d00,d01,d10,d11 =  wlinear_fit(x,y,w) ;
intercept = 0;
dintercept = 0.002;
slope = 6.77;
dslope = 0.25;

yfit = intercept + slope*x;
yup = intercept + dintercept + (slope+dslope)*x;
ydown = intercept - dintercept + (slope-dslope)*x;


plt.errorbar(x,y,yerr=deltax,fmt='r.',label='Data',color="blue")
plt.text(0.012,0.06822, r'Calculated h = $6.77 \times 10^{-34} Js $')
plt.text(0.012814,0.06414, r'$ \pm 0.25 \times 10^{-34} Js $')
plt.plot(x,yfit,label='Error weighted fit',color="orange")

plt.plot(x,yup,linestyle="dashed",linewidth=0.5,color="orange")
plt.plot(x,ydown,linestyle="dashed",linewidth=0.5,color="orange")

plt.legend()

plt.xlabel(r'$ \frac{c}{2deV}$ ( in $10^{-34}$)/(Js)');
plt.ylabel(r'$sin( \theta) $')


plt.grid()
plt.show();

