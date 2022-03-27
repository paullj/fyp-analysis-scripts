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
    df["Group"] = samples["Group"][i]
    df["Name"] = samples["Name"][i]

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

# %%
