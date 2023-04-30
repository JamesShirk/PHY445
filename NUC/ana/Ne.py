from scipy.integrate import quad
import numpy as np

def func(x, r, R):
    return 8*np.sqrt(r**2 - x**2)*np.sqrt(R**2 - x**2)
    return 8*r**3 * (np.cos(x)**2) * np.sqrt((R/r)**2 - np.sin(x)**2)

def ratio(r, R, h):
    volume = quad(func, 0, R, (r, R))
    #print(volume[0]/(np.pi*R**2*h))
    return volume[0]

def electrons(r = 1.5 * 2.54, R=3.25, h=10, rho = 2.7, Z=13, M = 1/(26.98)):
    V = ratio(r, R, h)
    return V*rho*M*6.022e23*Z*2

if __name__ == "__main__":
    print(electrons())