import numpy as np

# in m^2
def kleinNishina(theta,E):
    rad = theta*np.pi/180
    Eprime = E / (1+ (E/511)*(1-np.cos(rad)))
    g = Eprime/E
    r_e = 2.82e-15
    return g**2 * (g + 1/g + np.sin(theta*np.pi/180)) * (r_e**2 / 2)

# in m^2
def thompson(theta):
    r_e = 2.82e-15
    return  (r_e**2 / 2)*(1+np.cos(theta * np.pi/180)**2)

def fitXS(theta, a, E):
    rad = theta*np.pi/180
    Eprime = E / (1+ (E/511)*(1-np.cos(rad)))
    g = Eprime/E

    return g**2 * (g + 1/g + np.sin(theta*np.pi/180)) * a