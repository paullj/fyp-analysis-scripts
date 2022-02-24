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

# %%
input_file = 'sample_log.csv'  # A file that points to data files for each sample
output_dir = "/Users/paullj/Library/CloudStorage/OneDrive-UniversityofBath/FYP/Experiments/Analysis"

min_freq = 40
max_freq = 110e6

# %%
samples = pd.read_csv(input_file)
sample_data = [pd.read_csv(file, skiprows=(21), header=None, delim_whitespace=True, usecols=(0, 1), names=['Frequency', 'Real Impedance'], nrows=201, engine=("python"))
               for file in samples["Path"]]
groups = list(samples['Group'].unique())

print(f"Loaded {len(samples)} samples in {len(groups)} group(s) from `{input_file}`")

# %%
for i, df in enumerate(sample_data, start=0):
    df['Imaginary Impedance'] = pd.read_csv(samples["Path"][i], sep="\t", skiprows=(
        228), header=None, usecols=(1,), nrows=201, engine=("python"))

# Process and calcuate data
for i, df in enumerate(sample_data):
    k = (samples["Thickness"][i]*10**-3)/(samples["Area"][i]*10**-6)
    df["Permittivity"] = (-df["Imaginary Impedance"])/(df["Real Impedance"]
                                                       ** 2+df["Imaginary Impedance"]**2)*(k/(0.00000000000885*2*3.14*df["Frequency"]))
    df["AC Conductivity"] = (df["Real Impedance"])/(
        df["Real Impedance"]**2+df["Imaginary Impedance"]**2)*k
    df["Tan Delta"] = (-df["Real Impedance"] /
                       df["Imaginary Impedance"])
    df["Phase Angle"] = np.degrees(
        np.arctan(df["Imaginary Impedance"]/df["Real Impedance"]))
    df["Capacitance"] = (df["Permittivity"]
                         * (8.8541878128*10**-12)*(samples["Area"][i]*10**-6))/(samples["Thickness"][i]/1000)
    df["Impedance"] = np.sqrt(df["Real Impedance"]**2+df["Imaginary Impedance"]**2)
    df["Thickness"] = samples["Thickness"][i]
    df["Area"] = samples["Area"][i]

# %%
processed_data_dir = f'{output_dir}/Processed Data/'
if not os.path.exists(processed_data_dir):
    os.makedirs(processed_data_dir)


def get_valid_filename(s):
    return "".join(x for x in s.replace(' ', '_') if (x.isalnum() or x in "._-"))


for i, df in enumerate(sample_data):
    df.to_csv(
        f"{processed_data_dir}{get_valid_filename(samples.iloc[i]['Name'].lower())}_data.csv")

print(f"Read and processed {len(sample_data)} files")

# %% Setup figures
figures_dir = f'{output_dir}/Figures/'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)


def plot_against_frequency(samples, min_freq, max_freq, column, label, override_bounds=True):
    marker_list = ["o", "s", "v", "*", "D", "^", "d", ">", "D", "<", "P",
                   "H", "o", "s", "P", "v", "*", "D", "^", "d", ">", "D", "<", "H"]
    color_list = ['#e41a1c', '#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3', '#999999', '#dede00']
    for j, group in enumerate(groups):
        samples_in_group = samples.index[(
            samples['Group'] == group) & (samples['Include'] == 'X')]
        if len(samples_in_group) == 0:
            continue
        fig = plt.figure(j, figsize=(7, 5), dpi=300)

        for i in samples_in_group:
            _df = sample_data[i]
            if override_bounds and not np.isnan(samples.iloc[i]["Min Frequency"]):
                min_freq = samples.iloc[i]["Min Frequency"]
            if override_bounds and not np.isnan(samples.iloc[i]["Max Frequency"]):
                max_freq = samples.iloc[i]["Max Frequency"]

            df = _df[(_df["Frequency"] > min_freq) &
                     (_df["Frequency"] < max_freq)]

            plt.plot(df["Frequency"], df[column], linewidth=0.5,
                     marker=marker_list[i % len(marker_list)], markersize=5, color=color_list[i % len(color_list)], label=samples.iloc[i]["Name"], zorder=3)
            plt.xscale("log")
            # add the legend (will default to 'best' location)
            plt.legend(loc='best', fontsize=10, frameon=False)
            plt.xlabel(r"Frequency (Hz)", fontsize=12)
            plt.ylabel(f"{label}", fontsize=12)
            plt.grid(True, linestyle='dashed', linewidth='0.3',
                     color='grey', alpha=0.8, zorder=1)

        plt.tick_params(direction='in', which='major', length=3,
                        bottom=True, top=True, left=True, right=True)
        plt.tick_params(direction='in', which='minor', length=2,
                        bottom=True, top=True, left=True, right=True)

        fig_name = get_valid_filename(f'{group}_{column}').lower()
        plt.savefig(f'{figures_dir}{fig_name}.png')


# %% Plot Permittivity
plot_against_frequency(samples, min_freq, max_freq,
                       "Permittivity", "Permittivity")

# %% Plot Capacitance
plot_against_frequency(samples, min_freq, max_freq,
                       "Capacitance", "Capacitance (F)")

# %% Real Impedance
plot_against_frequency(samples, min_freq, max_freq,
                       "Real Impedance", "Real Impedance (Ω)")

# %% Imaginary Impedance
plot_against_frequency(samples, min_freq, max_freq,
                       "Imaginary Impedance", "Imaginary Impedance (Ω)", override_bounds=False)

# %% Plot AC Conductivity
plot_against_frequency(samples, min_freq, max_freq,
                       "AC Conductivity", "AC Conductivity (S/m)")

# %% Plot Phase Angle
plot_against_frequency(samples, 1e5, max_freq,
                       "Phase Angle", r"Phase angle ($^o$)", override_bounds=False)

# %% Plot Tanδ
plot_against_frequency(samples, min_freq, max_freq, "Tan Delta", "Tanδ")

# %% Plot Impedance
plot_against_frequency(samples, min_freq, max_freq,
                       "Impedance", "Impedance (Ω)", override_bounds=False)

# %%
