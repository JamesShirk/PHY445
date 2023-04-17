import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def eff():

    domega = (1.5**2 * np.pi) / ((11**2)*4*np.pi)

    Energy = np.array([10, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000])
    Efficiency = np.array([1, 1, 1, 1, 1, 0.99, 0.97, 0.96, 0.95, 0.93, 0.90, 0.87, 0.85, 0.83, 0.81])
    interp_eff = CubicSpline(Energy, Efficiency)

    energy_p2t = np.array([0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]) * 1000
    peak_to_total = np.array([0.95, 0.925, 0.8375, 0.725, 0.65, 0.575, 0.525, 0.4875, 0.45, 0.425])
    interp_p2t = CubicSpline(energy_p2t, peak_to_total)

    # x = np.logspace(1, 3, 1000)
    # plt.plot(energy_p2t, peak_to_total, "o")
    # plt.plot(x, interp_p2t(x), "--")
    # #plt.yscale('log')
    # #plt.xscale('log')
    # plt.xlim([150, 1000])
    # plt.grid()
    # plt.savefig("test.png")

    return interp_eff, interp_p2t