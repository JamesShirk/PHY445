import fitPeak as fp
import numpy as np
import matplotlib.pyplot as plt
# 'Pb','Fe','In','Brass','Nb','Au','V','Ag','Zn','Zr','Ti','Rh','Ni','W','U','Cu'

peaks = [
    [10.55, 12.61, 14.76], #Pb 10.28, 12.30, 14.30 La, Lb, Lg
    [6.40, 7.06], #Fe 6.76, 7.45 Ka, Kb
    [24.21, 27.28, 3.286], #In 
    #[8.05, 8.64, 8.91, 9.57], #Brass # to do
    [16.62, 18.62, 2.17, 2.26, 2.46], #Nb 16.63, 18.15 Ka, Kb
    [9.71, 11.44, 13.38], #Au 9.98, 11.94, 13.93 La, Lb
    [4.95, 5.43], #V 4.76, 5.22 Ka, Kb
    [22.16, 24.94, 2.98, 3.15, 3.52], #Ag 3.24 La
    [8.64, 9.57], #Zn 8.88, 9.85 Ka, Kb
    [15.78, 17.67], #Zr 16.02, 17.53 Ka, Kb
    [4.51, 4.93], #Ti 4.43, 4.84 Ka Kb
    [20.22, 22.72, 2.70, 2.83, 3.14], #Rh 
    [7.48, 8.26], #Ni 7.78, 8.61 Ka, Kb
    [8.40, 9.67, 11.29], #W # 8.35, 9.49, 10.97 La, Lb, Lg
    [13.61, 17.22, 20.17], #U 13.98, 17.10 La, Lb
    [8.05, 8.91], #Cu #8.16, 9.04 Ka, Kb
    [7.42, 8.40, 9.78] #yb #7.57, 8.59, 10.00 
]

# brass not included
AtomicNumbers = [82, 26, 49, 41, 79, 23, 47, 30, 40, 22, 45, 28, 74, 92, 29]


def moseley():
    # fe, nb, v, zn, zr, ti, ni, cu
    k_atomicNumbers = np.array([26, 41, 23, 30, 40, 22, 28, 29], dtype = np.float64)
    #k_lines = np.array([[6.40, 7.06], [16.62, 18.62], [4.95, 5.43], [8.64, 9.57], [15.78, 17.67], [4.51, 4.93], [7.48, 8.26], [8.05, 8.91]])
    k_lines_measured = np.array([[6.76, 7.45], [16.73, 18.15], [4.76, 5.22], [8.88, 9.85], [16.02, 17.53], [4.43, 4.84], [7.78, 8.61], [8.16, 9.04]])
    k_lines_error = np.array([[0.33, 0.43], [0.17, 0.21], [0.08, 0.11], [0.21, 0.24], [0.29, 0.16], [0.11, 0.12], [0.23, 0.28], [0.16, 0.17]])
    # pb, au, ag, w, u
    #l_atomicNumbers = np.array([82, 79, 47, 74, 92])
    l_atomicNumbers = np.array([82, 79, 74, 92, 70, 47])
    l_lines = np.array([[10.55, 12.61, 14.76], [9.71, 11.44], [8.40, 9.67, 11.29], [13.61, 17.22], [7.42, 8.40, 9.78]], dtype = 'object')
    # ag = 47, 3.24, 1.08
    l_lines_measured = np.array([[10.28, 12.30, 14.30], [9.99, 11.94, 13.93], [8.91, 9.49, 10.97], [13.98, 17.10], [7.57, 8.59, 10.00 ], [3.24]], dtype = 'object')
    l_lines_error = np.array([[0.22, 0.1, 0.19], [0.41, 0.14, 1.27],[0.58, 0.34, 1.08], [1.30, 5.77], [0.08, 0.25, 0.27], [1.08]], dtype = 'object')
    # read in values from arrays above
    K_alpha = np.empty_like(k_atomicNumbers)
    K_alpha_error = np.empty_like(k_atomicNumbers)
    K_alpha_sys = np.empty_like(k_atomicNumbers)
    K_beta = np.empty_like(k_atomicNumbers)
    K_beta_error = np.empty_like(k_atomicNumbers)
    K_beta_sys = np.empty_like(k_atomicNumbers)

    a = 4.88
    vara = 0.85
    b = 34.29
    varb = 8.66



    for i in range(len(k_lines_measured)):
        K_alpha[i] = k_lines_measured[i][0]
        K_alpha_chan = k_lines_measured[i][0]*b + a
        K_alpha_sys[i] = np.abs((1/b)*((k_lines_error[i][0] - vara) - (varb/b)*(K_alpha_chan - a)))
        K_alpha_error[i] = k_lines_error[i][0]
        K_beta[i] = k_lines_measured[i][1]
        K_beta_error[i] = k_lines_error[i][1]
        K_beta_chan = k_lines_measured[i][1]*b + a
        K_beta_sys[i] = np.abs((1/b)*((k_lines_error[i][1] - vara) - (varb/b)*(K_beta_chan - a)))

    #sort
    idx = np.argsort(k_atomicNumbers)
    K_alpha = K_alpha[idx]
    K_alpha_error = K_alpha_error[idx]
    K_alpha_error = np.divide(K_alpha_error, 2*np.sqrt(K_alpha))
    K_beta = K_beta[idx]
    K_beta_error = K_beta_error[idx]
    K_beta_error = np.divide(K_beta_error, 2*np.sqrt(K_beta))
    k_atomicNumbers = k_atomicNumbers[idx]

    L_lines = [[], [], []]
    L_error = [[], [], []]
    L_Z = [[], [], []]
    L_sys = [[], [], []]

    for i in range(len(l_lines_measured)):
        if len(l_lines_measured[i]) >= 1:
            L_lines[0].append(l_lines_measured[i][0])
            L_error[0].append(l_lines_error[i][0])
            L_chan = l_lines_measured[i][0]*b + a
            L_sys[0].append(np.abs((1/b)*((l_lines_error[i][0] - vara) - (varb/b)*(L_chan - a))))
            L_Z[0].append(l_atomicNumbers[i])
        if len(l_lines_measured[i]) >= 2:
            L_lines[1].append(l_lines_measured[i][1])
            L_error[1].append(l_lines_error[i][1])
            L_chan = l_lines_measured[i][0]*b + a
            L_sys[1].append(np.abs((1/b)*((l_lines_error[i][1] - vara) - (varb/b)*(L_chan - a))))
            L_Z[1].append(l_atomicNumbers[i])
        if len(l_lines_measured[i]) == 3:
            L_lines[2].append(l_lines_measured[i][2])
            L_error[2].append(l_lines_error[i][2])
            L_chan = l_lines_measured[i][0]*b + a
            L_sys[2].append(np.abs((1/b)*((l_lines_error[i][2] - vara) - (varb/b)*(L_chan - a))))
            L_Z[2].append(l_atomicNumbers[i])

    

    L_lines = np.array(L_lines, dtype = 'object')
    L_error = np.array(L_error, dtype = 'object')
    L_sys = np.array(L_sys, dtype = 'object')
    for i in range(len(L_error)):
        idx = np.argsort(L_Z[i])
        L_Z[i] = np.array(L_Z[i])[idx]
        L_lines[i] = np.array(L_lines[i])[idx]
        L_error[i] = np.array(L_error[i])[idx]
        L_sys[i] = np.array(L_sys[i])[idx]
        L_error[i] = L_error[i]/(2*np.sqrt(L_lines[i]))
    L_Z = np.array(L_Z, dtype = 'object')
    # weightd linear fits
    fit = np.empty((5, 6), dtype = np.float64)
    fit[0] = fp.wlinear_fit(k_atomicNumbers, np.sqrt(K_alpha), 1/(K_alpha_error**2))
    fit[1] = fp.wlinear_fit(k_atomicNumbers, np.sqrt(K_beta), 1/(K_beta_error**2))
    fit[2] = fp.wlinear_fit(L_Z[0], np.sqrt(L_lines[0]), 1/(L_error[0]**2))
    fit[3] = fp.wlinear_fit(L_Z[1], np.sqrt(L_lines[1]), 1/(L_error[1]**2))
    fit[4] = fp.wlinear_fit(L_Z[2], np.sqrt(L_lines[2]), 1/(L_error[2]**2))
    xfit_K = np.linspace(k_atomicNumbers[0], k_atomicNumbers[-1], 100)
    xfit_L = np.empty((3, 100))
    for i in range(len(xfit_L)):
        xfit_L[i] = np.linspace(L_Z[i][0], L_Z[i][-1], 100)
    # plotting
    fig, axs = plt.subplots(2, 1, figsize = (10, 10))
    axs[0].errorbar(k_atomicNumbers, np.sqrt(K_alpha), yerr = K_alpha_error, marker = "o", ls = 'None', label = r"$K_{\alpha}$")
    axs[0].plot(xfit_K, fit[0][1]*xfit_K + fit[0][0], "--", label = r"$K_{\alpha}$ fit")#, $\chi^2$ = {}".format(fit[0][-1]))
    axs[0].errorbar(k_atomicNumbers, np.sqrt(K_beta), yerr = K_beta_error, marker = "o", ls = 'None', label = r"$K_{\beta}$")
    axs[0].plot(xfit_K, fit[1][1]*xfit_K + fit[1][0], "--", label = r"$K_{\beta}$ fit")#, $\chi^2$ = {}".format(fit[1][-1]))
    axs[0].grid()
    axs[0].legend()
    axs[0].text(35.2, 2.8, r"$\chi^2_{\alpha}\ /\ \mathrm{DoF}$" + " = {0}".format(round(fit[0][-1]/7, 2)))
    axs[0].text(35.2, 2.6,r"$\chi^2_{\beta}\ /\ \mathrm{DoF}$" + " = {0}".format(round(fit[1][-1]/7, 2)))
    
    lab = [r"$L_{\alpha}$", r"$L_{\beta}$", r"$L_{\gamma}$"]
    for i in range(len(L_Z)):
        axs[1].errorbar(L_Z[i], np.sqrt(L_lines[i]), yerr = L_error[i], marker = "o", ls = 'None', label = lab[i])
        axs[1].plot(xfit_L[i], fit[i+2][1]*xfit_L[i] + fit[i+2][0], "--", label = lab[i] + r" fit")

    axs[1].legend()
    axs[1].grid()
    axs[1].text(81, 2.5, r"$\chi^2_{\alpha}\ /\ \mathrm{DoF}$" + " = {0}".format(round(fit[2][-1]/5, 2)))
    axs[1].text(81, 2.3,r"$\chi^2_{\beta}\ /\ \mathrm{DoF}$" + " = {0}".format(round(fit[3][-1]/4, 2)))
    axs[1].text(81, 2.1,r"$\chi^2_{\gamma}\ /\ \mathrm{DoF}$" + " = {0}".format(round(fit[4][-1]/3, 2)))
    fig.suptitle("Moseley's law")
    fig.supxlabel("Z")
    fig.supylabel(r"$\sqrt{E \ (keV)}$")
    plt.tight_layout()
    plt.savefig("plots/finalMoseley.png")
    plt.close()

    nam = ["ka", "kb", "la", "lb", "lg"]
    a = [3/4, 8/9, 5/36, 5/36, 3/16]
    for i in range(len(nam)):
        print("Screening constant for " + nam[i] + " is {0} +- {1}".format(
            round(-1*fit[i][0]*np.sqrt(1000/(13.6*a[i])), 2),
            round(np.sqrt(fit[i][2])*np.sqrt(1000/(13.6*a[i])), 2)
        ))


    fit_sys = np.empty((5, 6), dtype = np.float64)
    fit_sys[0] = fp.wlinear_fit(k_atomicNumbers, np.sqrt(K_alpha), 1/(K_alpha_sys**2))
    fit_sys[1] = fp.wlinear_fit(k_atomicNumbers, np.sqrt(K_beta), 1/(K_beta_sys**2))
    fit_sys[2] = fp.wlinear_fit(L_Z[0], np.sqrt(L_lines[0]), 1/(L_sys[0]**2))
    fit_sys[3] = fp.wlinear_fit(L_Z[1], np.sqrt(L_lines[1]), 1/(L_sys[1]**2))
    fit_sys[4] = fp.wlinear_fit(L_Z[2], np.sqrt(L_lines[2]), 1/(L_sys[2]**2))
    print()

    nam = ["ka", "kb", "la", "lb", "lg"]
    a = [3/4, 8/9, 5/36, 5/36, 3/16]
    for i in range(len(nam)):
        print("Screening constant for " + nam[i] + " is {0} +- {1}".format(
            round(-1*fit_sys[i][0]*np.sqrt(1000/(13.6*a[i])), 2),
            round(np.sqrt(fit_sys[i][2])*np.sqrt(1000/(13.6*a[i])), 2)
        ))

if __name__ == "__main__":
    moseley()