# %% Import libraries
import numpy as np
import matplotlib.pyplot as plt

# %% Define variables
# Sweep settings
min_freq = 1e7
max_freq = 110e7
num_points = 201

# Material Properties
eS_33r = -3.36391 # Relative dielectric constant
e_0 = 8.854e-12 # Dielectric constant of vacuum
C_33 = 194e9 # Elastic stiffness coefficient
E_33 = 5.11e9 # Piezoelectric coefficient 
rho = 5700 # Density

delta_m = 0.1 # Mechanical losses
delta_e = 0.1 # Dielectric losses

# Material Geometry
A = 97.022e-6 # Electrode area
t = 0.047e-3 # Thickness

# %% Calculate impedance
k_t = E_33 / np.sqrt(C_33*eS_33r, dtype=complex)
v_t = np.sqrt(C_33 / rho, dtype=complex)

delta = (1-k_t**2)*delta_m + (k_t**2)*delta_e

C_33s = C_33 * (1+ 1j * delta_m)
eS_33s = C_33 * (1+ 1j * delta_e)

k_ts2 = (k_t**2)/((1+1j*delta)*(1-1j*delta_e))
v_ts = v_t * np.sqrt(1+1j*delta)

freq = np.linspace(min_freq, max_freq, num_points)
omega = 2 * np.pi * freq

Ys = ((1j*omega*eS_33r*e_0*A)/t) * ((1-k_ts2) * 1/((np.tan(omega*t/(2*v_ts))/(omega*t/(2*v_ts)))))
Zs = 1 / Ys

# %% Plot results
fig = plt.figure()
plt.plot(freq[1:], Zs.real[1:], label="Real Component")
plt.plot(freq[1:], Zs.imag[1:], label="Imaginary Component")
plt.xscale("log")

plt.legend(loc='best', fontsize=10, frameon=False)
plt.xlabel(r"Frequency (Hz)", fontsize=12)
plt.ylabel(f"Impedance (Ω)", fontsize=12)
plt.grid(True, linestyle='dashed', linewidth='0.3',
            color='grey', alpha=0.8, zorder=1)

plt.tick_params(direction='in', which='major', length=3,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)

# %%
