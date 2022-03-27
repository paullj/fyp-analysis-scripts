# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 2022

@author: Paul Lavender-Jones
"""
# %% Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# %%

csv_to_read = 'candidate_substrates.csv'
materials = pd.read_csv(csv_to_read)

t = np.linspace(0.1e-3, 3e-3, 55)
R = 0
k_iq = 0.8
n = 1

fig = plt.figure(figsize=(6,3), dpi=300)

for i, name in enumerate(materials["Name"]):
    c_33 = materials["Stiffness"][i]
    rho = materials["Density"][i]
    f = 1/(2*t)*(1-R-(4*k_iq**2)/np.pi**2)*np.sqrt(c_33**n/rho)
    plt.plot(t * 1e3, f * 1e-6, label=name)

plt.xscale("log")
plt.ylabel("Resonant Frequency (MHz)")
plt.xlabel("Disc Thickness (mm)")
plt.legend(loc='upper right', fontsize=10, frameon=True)
plt.grid(True, linestyle='dashed', linewidth='0.3',
                color='grey', alpha=0.8, zorder=1)
plt.tick_params(direction='in', which='major', length=3,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)

# %%
