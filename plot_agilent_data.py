# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 2022 at 12:34pm 

@author: Paul Lavender-Jones
"""
# %% Import libraries
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import get_material_impedance as model

%matplotlib widget

# %% 
min_freq = 10e3
max_freq = 110e6 

path = '/Users/paullj/Library/CloudStorage/OneDrive-UniversityofBath/FYP/Experiments/Analysis/Processed Data/ref_data.csv'
output_dir = "/Users/paullj/Library/CloudStorage/OneDrive-UniversityofBath/FYP/Experiments/Analysis"

# %% 
df = pd.read_csv(path)

df = df[(df["Frequency"] > min_freq) &
        (df["Frequency"] < max_freq)]

def get_valid_filename(s):
    return "".join(x for x in s.replace(' ', '_') if (x.isalnum() or x in "._-"))

figures_dir = f'{output_dir}/Figures/'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

name = df.iloc[0]["Name"]

# %% Plot |Z| and Phase Angle
fig = plt.figure()

z_plot = plt.subplot(211)
plt.plot(df["Frequency"], df["Impedance"], label=name, color='#e41a1c', zorder=3)
plt.yscale("log")
plt.xscale("log")
plt.legend(loc='lower left', fontsize=10, frameon=True)
plt.ylabel(f"|Z| (Ω)", fontsize=12)
plt.grid(True, linestyle='dashed', linewidth='0.3',
            color='grey', alpha=0.8, zorder=1)
plt.tick_params(direction='in', which='major', length=3,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)

theta_plot = plt.subplot(212)

plt.plot(df["Frequency"], df["Phase Angle"],label=name, color='#377eb8', zorder=3)
plt.xscale("log")
plt.legend(loc='upper left', fontsize=10, frameon=True)
plt.xlabel(r"Frequency (Hz)", fontsize=12)
plt.ylabel(r"$\Theta$ ($^o$)", fontsize=12)
plt.grid(True, linestyle='dashed', linewidth='0.3',
            color='grey', alpha=0.8, zorder=1)
plt.tick_params(direction='in', which='major', length=3,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)

fig_name = get_valid_filename(f'{name}_z_theta').lower()
plt.savefig(f'{figures_dir}{fig_name}.png')

# %% Add expected resonance values
c_33 = 60e9 # Elastic stiffness coefficient (N/m^2)
rho = 7500 # Density (kg/m^3)
A = df.iloc[0]["Area"] * 1e-6          # Electrode area (m^2)
t = df.iloc[0]["Thickness"] * 1e-3        # Thickness (m)

R = 0
k_iq = 0
for n in range(1, 2):
    f = 1/(2*t)*(1-R-(4*k_iq**2)/np.pi**2)*np.sqrt(c_33**n/rho)
    z_plot.axvline(x=f, linestyle="-.", color="gray", zorder=0)

# %% Add simulated impedance to plot
num_points = 20001    # Number of points between start and end

epsilon_33 = 1200 #173 #400    # Relative dielectric constant
d_33 = 515               # Piezoelectric coefficient (pC/N)

delta_m = 2            # Mechanical losses
delta_e = 0.08           # Dielectric losses

fixture_weight = 46.02*9.81
fixture_stress = fixture_weight/ (np.pi*2**2)
e_33 = d_33 * fixture_stress

(freq, Z, f_1) = model.get_impedance(epsilon_33, c_33, e_33, rho, A, t, delta_e, delta_m, min_freq=min_freq, max_freq=max_freq, num_points=num_points)

z_plot.plot(freq[1:], np.sqrt(Z.real[1:]**2+Z.imag[1:]**2), label="Model", linestyle="--", color='#ff8c8e')
z_plot.legend(loc='lower left', fontsize=10, frameon=True)
theta_plot.plot(freq[1:], np.degrees(np.arctan(Z.imag[1:]/Z.real[1:])), label="Model", linestyle="--", color='#8ce5ff')
theta_plot.legend(loc='upper left', fontsize=10, frameon=True)
                        
# %%
fig = plt.figure()
plt.plot(df["Frequency"], df["Real Impedance"], zorder=3)
plt.plot(freq[1:], Z.real[1:])
plt.xscale("log")

# %%
plt.savefig(f'{figures_dir}{fig_name}.png')

# %%
