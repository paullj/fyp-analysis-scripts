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

input_file = 'samples_to_model.csv'  # A file that points to data files for each sample

min_freq = 10e3 # 10 kHz
max_freq = 110e6 # 110 MHz

# %%
samples = pd.read_csv(input_file)
sample_data = [pd.read_csv(file) for file in samples["Path"]]

# %%
for i, df in enumerate(sample_data):
    df = df[(df["Frequency"] > min_freq) &
        (df["Frequency"] < max_freq)]

    fig = plt.figure()

    z_plot = plt.subplot(211)
    plt.plot(df["Frequency"], df["Impedance"], label=samples.iloc[i]["Name"], color='#e41a1c', zorder=3)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(min_freq, max_freq)

    plt.legend(loc='lower left', fontsize=10, frameon=True)
    plt.ylabel(f"|Z| (Ω)", fontsize=12)
    plt.grid(True, linestyle='dashed', linewidth='0.3',
                color='grey', alpha=0.8, zorder=1)
    plt.tick_params(direction='in', which='major', length=3,
                    bottom=True, top=True, left=True, right=True)
    plt.tick_params(direction='in', which='minor', length=2,
                            bottom=True, top=True, left=True, right=True)

    theta_plot = plt.subplot(212)

    plt.plot(df["Frequency"], df["Phase Angle"],label=samples.iloc[i]["Name"], color='#377eb8', zorder=3)

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
    plt.xlim(min_freq, max_freq)

    c_33 = samples.iloc[i]["Stiffness"]
    rho = samples.iloc[i]["Density"]
    A = samples.iloc[i]["Area"]
    t = samples.iloc[i]["Thickness"]

    R = 0
    k_iq = 0
    for n in range(1, 4):
        f = 1/(2*t)*(1-R-(4*k_iq**2)/np.pi**2)*np.sqrt(c_33**n/rho)
        z_plot.axvline(x=f, linestyle="-.", color="gray", zorder=0)
        theta_plot.axvline(x=f, linestyle="-.", color="gray", zorder=0)

    num_points = 20001    # Number of points between start and end

    epsilon_33 = 1200 #173 #400    # Relative dielectric constant
    d_33 = 515               # Piezoelectric coefficient (pC/N)

    delta_m = 0.50            # Mechanical losses
    delta_e = 0.05           # Dielectric losses

    fixture_weight = 46.02*9.81
    fixture_stress = fixture_weight/ (np.pi*2**2)
    e_33 = d_33 * fixture_stress

    (freq, Z, f_1) = model.get_impedance(epsilon_33, c_33, e_33, rho, A, t, delta_e, delta_m, min_freq=min_freq, max_freq=max_freq, num_points=num_points)

    z_plot.plot(freq[1:], np.sqrt(Z.real[1:]**2+Z.imag[1:]**2), label="Model", linestyle="--", color='#ff8c8e')
    z_plot.legend(loc='lower left', fontsize=10, frameon=True)
    theta_plot.plot(freq[1:], np.degrees(np.arctan(Z.imag[1:]/Z.real[1:])), label="Model", linestyle="--", color='#8ce5ff')
    theta_plot.legend(loc='upper left', fontsize=10, frameon=True)

#     BCZT Ceramic,1.88,76.553,7500,60000000000
# BT Polymer,0.043,182.605,7500,60000000000
# %%
