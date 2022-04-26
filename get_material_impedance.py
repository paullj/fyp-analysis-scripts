import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt

# Using equation in 'Comparison of several methods to characterise the high frequency behaviour of piezoelectric ceramics for transducer applications'
# From here: https://www.sciencedirect.com/science/article/pii/S0041624X99000591

def get_impedance(eps_33, c_33, e_33, rho, A, t, delta_m, delta_e, min_freq=40, max_freq=110e7, num_points=201):
    e_0 = 8.8541878e-12 # Dielectric constant of vacuum (F/m)

    # %% Calculate impedance
    k_t = e_33 / np.sqrt(c_33*eps_33)
    v_t = np.sqrt(c_33 / rho)

    delta = (1-k_t**2)*delta_m + (k_t**2)*delta_e

    c_33s = c_33 * (1 + 1j*delta_m)
    eps_33s = eps_33 * (1 - 1j*delta_e)

    k_ts2 = (k_t**2)/((1+1j*delta)*(1-1j*delta_e))
    v_ts = v_t * np.sqrt(1+1j*delta)

    freq = np.linspace(min_freq, max_freq, num_points)
    omega = 2 * np.pi * freq

    Y = ((1j*omega*eps_33s*e_0*A)/t) * ((1-k_ts2)*((np.tan(omega*t/(2*v_ts))/(omega*t/(2*v_ts)))))**-1
    Z = 1 / Y

    # fc = lambda omega: 1- k_t**2 / ((omega*t/(2*v_t))*(1/(np.tan(omega*t/(2*v_t)) - 0*omega*t/(2*v_t))))
    # x = opt.fsolve(fc, np.array([1e5, 3e5, 6e5]))
    # print(x)

    R = 0
    f_1 = 1/(2*t)*(1-R-(4*k_t**2/(np.pi**2)))*(c_33/rho)**0.5
    
    return (freq, Z, f_1)