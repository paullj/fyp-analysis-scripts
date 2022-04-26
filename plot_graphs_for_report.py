# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 2022 at 12:34pm 

@author: Paul Lavender-Jones
"""
# %% Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-colorblind')

output_dir = "/Users/paullj/Library/CloudStorage/OneDrive-UniversityofBath/FYP/Documents/Final Report/Figures"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#%% Plot survival rate vs treatment cost
fig = plt.figure(figsize=(7,5), dpi=300)

df = pd.DataFrame({
 'survival_rate' : [81.11,76.11,59.99,29.99,17.22,7.77],
 'cost' : [0,26.24,64.17,105.97,125.42,140]
 })

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

left = df['survival_rate'].plot(kind='bar')
left.set_ylabel(f"Survival Rate (%)", fontsize=12)
left.set_ylim(0, 100)
left.set_xticklabels(["No Cancer","Stage 0","Stage I","Stage II","Stage III","Stage IV"])
left.set_xlabel(f"Cancer Stage at Diagnosis", fontsize=12)

left.tick_params(direction='in', which='major', length=3,
                bottom=False, top=False, left=True, right=False)
left.tick_params(direction='in', which='minor', length=2,
                        bottom=False, top=False, left=True, right=False)

right = df['cost'].plot(color=colors[1], marker='o', secondary_y=True)

right.set_ylabel(f"Cost (Thousands)", fontsize=12)

right.yaxis.set_major_formatter('${x:1.0f}')
right.tick_params(direction='in', which='major', length=3,
                bottom=False, top=False, left=False, right=True)
right.tick_params(direction='in', which='minor', length=2,
                        bottom=False, top=False, left=False, right=True)

h1, l1 = left.get_legend_handles_labels()
h2, l2 = right.get_legend_handles_labels()

lgnd = plt.legend(h1+h2,["Proportion alive 15 years after diagnosis", "Total treatment cost"], loc='upper left', fontsize=10, frameon=True)
# plt.savefig(f'{output_dir}/survival_rate_treatment_cost.png')
plt.show()

#%% Plot survival rate vs treatment cost
fig = plt.figure(figsize=(7,5), dpi=300)
survival_rate = [100, 100, 100, 95.6, 49.0]

plt.bar(["Stage 0","Stage I","Stage II","Stage III","Stage IV"],survival_rate)
plt.ylabel(f"Survival Rate (%)", fontsize=12)
# left.set_xticklabels(["No Cancer","Stage 0","Stage I","Stage II","Stage III","Stage IV"])
plt.xlabel(f"Cancer Stage at Diagnosis", fontsize=12)

plt.tick_params(direction='in', which='major', length=3,
                bottom=False, top=False, left=True, right=False)
plt.tick_params(direction='in', which='minor', length=2,
                        bottom=False, top=False, left=True, right=False)

# right = df['cost'].plot(color=colors[1], marker='o', secondary_y=True)

# right.set_ylabel(f"Cost (Thousands)", fontsize=12)

# right.yaxis.set_major_formatter('${x:1.0f}')
# right.tick_params(direction='in', which='major', length=3,
#                 bottom=False, top=False, left=False, right=True)
# right.tick_params(direction='in', which='minor', length=2,
#                         bottom=False, top=False, left=False, right=True)

# h1, l1 = left.get_legend_handles_labels()
# h2, l2 = right.get_legend_handles_labels()
plt.savefig(f'{output_dir}/survival_rate_treatment_cost.png')
plt.show()

# %%
total_diagnoses = [
    143602,146113,151872,156723,157383,158302,161456,161313,164994,170282,174943,176109,180209,179985,190322,188066,192850,197703,198242,199299,202715,211169,208096,213734,214998,215528,222288,223387,230827,230727,234558,232756,237072,243067,246684,253480,255288,266590,272622,275908,281472,289656,298922,299983,299949,303388,306000
]
years = [
    1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017
]
fig = plt.figure(figsize=(7,5), dpi=300)

ax = plt.plot(years, total_diagnoses, label="Amount of people diagnosed with cancer",marker="o")
current_values = plt.gca().get_yticks()
plt.gca().set_yticklabels(['{:,.0f}'.format(x) for x in current_values])

plt.legend(loc='best', fontsize=10, frameon=True)
plt.xlabel("Year", fontsize=12)
plt.ylabel("Total Diagnoses", fontsize=12)
plt.grid(True, linestyle='dashed', linewidth='0.3',
    color='grey', alpha=0.8, zorder=1)

plt.tick_params(direction='in', which='major', length=3,
    bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in', which='minor', length=2,
    bottom=True, top=True, left=True, right=True)
plt.savefig(f'{output_dir}/growing_cancer_diagnoses.png')


# %%
csv_to_read = 'candidate_substrates.csv'
materials = pd.read_csv(csv_to_read)

t = np.geomspace(0.1e-3, 3e-3, 20)
R = 0
k_iq = 0.8
n = 1

fig = plt.figure(figsize=(6,3), dpi=300)

marker_list = ["o", "s", "v", "*", "D", "^", "d", ">", "D", "<", "P",
                   "H", "o", "s", "P", "v", "*", "D", "^", "d", ">", "D", "<", "H"]

for i, name in enumerate(materials["Name"]):
    c_33 = materials["Stiffness"][i]
    rho = materials["Density"][i]
    # f = 1/(2*t)*(1-R-(4*k_iq**2)/np.pi**2)*np.sqrt(c_33**n/rho)
    # f = 1/(2*t)*(1-(4*k_iq**2)/np.pi**2)*np.sqrt(c_33/rho)
    f = np.sqrt(c_33/rho)/(2*t)
    plt.plot(t * 1e3, f * 1e-6, label=name, marker=marker_list[i])

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
plt.savefig(f'{output_dir}/thickness_res_freq.png')

# %%
