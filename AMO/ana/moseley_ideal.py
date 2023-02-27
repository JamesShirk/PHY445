import numpy as np
import fitPeak as fp
from scipy.stats import linregress

k_atomicNumbers = np.array([26, 41, 23, 30, 40, 22, 28, 29], dtype = np.float64)
k_lines = np.array([[6.40, 7.06], [16.62, 18.62], [4.95, 5.43], [8.64, 9.57], [15.78, 17.67], [4.51, 4.93], [7.48, 8.26], [8.05, 8.91]])
idx = np.argsort(k_atomicNumbers)
k_alpha = k_lines.T[0][idx]
k_beta = k_lines.T[1][idx]
k_atomicNumbers = k_atomicNumbers[idx]
print("kalpha")
print(-1*linregress(k_atomicNumbers, np.sqrt(k_alpha))[1]*np.sqrt(1000/(13.6*(3/4))))
print("kbeta")
print(-1*linregress(k_atomicNumbers, np.sqrt(k_beta))[1]*np.sqrt(1000/(13.6*(8/9))))

l_lines = np.array([[10.55, 12.61, 14.76], [9.71, 11.44], [8.40, 9.67, 11.29], [13.61, 17.22], [7.42, 8.40, 9.78], [2.98]], dtype = 'object')
l_alpha = np.array([10.55, 9.71, 8.40, 13.61, 7.42, 2.98])
l_alpha_z = np.array([82, 79, 74, 92, 70, 47])
l_beta = np.array([12.61, 11.44, 9.67, 17.22, 8.40])
l_beta_z = np.array([82, 79, 74, 92, 70])
l_gamma = np.array([14.76, 13.38, 11.29, 9.78])
l_gamma_z = np.array([82, 79, 74, 70])
print("lalpha")
print(-1*linregress(l_alpha_z, np.sqrt(l_alpha))[1]*np.sqrt(1000/(13.6*(5/36))))
print("lbeta")
print(-1*linregress(l_beta_z, np.sqrt(l_beta))[1]*np.sqrt(1000/(13.6*(5/36))))
print("lgamma")
print(-1*linregress(l_gamma_z, np.sqrt(l_gamma))[1]*np.sqrt(1000/(13.6*(3/16))))