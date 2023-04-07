#calculate t2star and t2'' 
import numpy as np
import matplotlib.pyplot as plt

names = ["cs","wt", "mo", "gy", "po", "rb"]
x = [1, 2, 3, 4, 5, 6]
t2p = np.array([20.7,907, 15.4, 31.9, 24.6, 19.3])
errort2p = np.array([0.2, 38, 0.3, 0.4, 0.5, 1.1])
t1 = np.array([20.2, 2842, 26.6, 39.6, 38.5, 26.6])
errort1 = np.array([0.1, 32, 0.3, 0.6, 0.5, 0.2])
t2 = np.array([0.62, 0.63, 0.18, 0.34, 0.31, 1.72])
errort2 = np.array([0.06, 0.05, 0.01, 0.01, 0.01, 0.17])
chiSuc = np.array([-9.035e-7])


T2_star = np.divide(np.multiply(t2p, t2), t2p-t2)
errorT2_star = T2_star**2 * np.sqrt(np.divide(errort2p, t2p**2)**2 + np.divide(errort2, t2**2)**2)
T2_2p = np.divide(2*np.multiply(t1, t2p), 2*t1-t2p)
error_T22p = T2_2p**2 * np.sqrt(np.divide(errort1, t1**2)**2 + 4*np.divide(errort2p, t2p**2)**2)

print(T2_star)
print(errorT2_star)
print(T2_2p)
print(error_T22p)

weighted_average = np.sum(np.dot(T2_star[0:-1],1/errorT2_star[0:-1]**2)) / (np.sum(1/errorT2_star[0:-1]**2))
weighted_average_error = np.sqrt(5/4)*np.sqrt(np.sum(np.dot(1/errorT2_star[0:-1]**2, (T2_star[0:-1]-weighted_average)**2)))/ (np.sum(1/errorT2_star[0:-1]**2))
print("Average = {0} pm {1}".format(
    round(weighted_average, 3), 
    round(weighted_average_error,3)
    ))

plt.errorbar(x, T2_star, yerr = errorT2_star, marker = "o", ls = 'None', fillstyle = 'none', label = "data")
plt.xticks(x, names, rotation = 'vertical')
plt.text(2, 1.5, r"Average excluding rubber = {0} $\pm$ {1} (ms)".format(
    round(weighted_average, 2), 
    0.01
    ))
plt.xlabel("Name")
plt.ylabel(r"$T_2^*$ (ms)")
plt.savefig("T2_star.png", dpi = 300)