# %% Import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import get_material_impedance as model

%matplotlib widget

# %% Define variables
# Sweep settings
min_freq = 1e3       # Sweep start frequency (Hz)
max_freq = 110e6       # Sweep end frequency (Hz)
num_points = 20001    # Number of points between start and end

# Material Properties
epsilon_33 = 4000        # Relative dielectric constant
c_33 = 9.85e10           # Elastic stiffness coefficient (N/m^2)
d_33 = 280               # Piezoelectric coefficient (pC/N)
rho = 5570               # Density (kg/m^3)

delta_m = 0.02           # Mechanical losses
delta_e = 0.08           # Dielectric losses

# Material Geometry
A = 78.719e-6          # Electrode area (m^2)
t = 1.84e-3            # Thickness (m)

fixture_weight = 46.02*9.81
fixture_stress = fixture_weight/ (np.pi*2**2)
e_33 = d_33 * fixture_stress
(freq, Z, f_1) = model.get_impedance(epsilon_33, c_33, e_33, rho, A, t, delta_e, delta_m, min_freq=min_freq, max_freq=max_freq, num_points=num_points)

path = '/Users/paullj/Library/CloudStorage/OneDrive-UniversityofBath/FYP/Experiments/Analysis/Processed Data/bczt_2_data.csv'
df = pd.read_csv(path)

df = df[(df["Frequency"] > min_freq) &
        (df["Frequency"] < max_freq)]

# %% Plot results
fig = plt.figure()
plt.axvline(x=f_1, linestyle="-.", color="gray")
plt.plot(freq[1:], Z.real[1:], label="Model Real", linestyle="--", color='#ff8c8e')
plt.plot(df["Frequency"], df["Real Impedance"], label="Real Component", color='#e41a1c', zorder=3)
plt.xscale("log")
plt.legend(loc='lower right', fontsize=10, frameon=False)
plt.xlabel(r"Frequency (Hz)", fontsize=12)
plt.ylabel(f"Impedance (Ω)", fontsize=12)
plt.grid(True, linestyle='dashed', linewidth='0.3',
            color='grey', alpha=0.8, zorder=1)
plt.tick_params(direction='in', which='major', length=3,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)

# %% 
fig = plt.figure()
plt.axvline(x=f_1, linestyle="-.", color="gray", zorder=0)
plt.plot(freq[1:], Z.imag[1:], label="Model Imaginary", linestyle="--", color='#8ce5ff')
plt.plot(df["Frequency"], df["Imaginary Impedance"], label="Imaginary Component", color='#377eb8', zorder=3)
plt.xscale("log")
plt.legend(loc='lower right', fontsize=10, frameon=False)
plt.xlabel(r"Frequency (Hz)", fontsize=12)
plt.ylabel(f"Impedance (Ω)", fontsize=12)
plt.grid(True, linestyle='dashed', linewidth='0.3',
            color='grey', alpha=0.8, zorder=1)
plt.tick_params(direction='in', which='major', length=3,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)
