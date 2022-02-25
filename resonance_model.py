# %% Import libraries
import numpy as np
import matplotlib.pyplot as plt
import get_material_impedance as model

# %% Define variables
# Sweep settings
# min_freq = 40       # Sweep start frequency (Hz)
# max_freq = 110e7    # Sweep end frequency (Hz)
num_points = 20001    # Number of points between start and end

# Material Properties
epsilon_33 = 3      # Relative dielectric constant
c_33 = 194e9        # Elastic stiffness coefficient (N/m^2)
e_33 = 5.11e12      # Piezoelectric coefficient (pC/m^2)
rho = 2928          # Density (kg/m^3)

delta_m = 0.10      # Mechanical losses
delta_e = 0.10      # Dielectric losses

# Material Geometry
A = 100e-6       # Electrode area (m^2)
t = 2e-3            # Thickness (m)

(freq, Z) = model.get_impedance(epsilon_33, c_33, e_33, rho, A, t, delta_e, delta_m, num_points=num_points)

# %% Plot results
fig = plt.figure()
plt.plot(freq[1:], Z.real[1:], label="Real Component")
plt.plot(freq[1:], Z.imag[1:], label="Imaginary Component")
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