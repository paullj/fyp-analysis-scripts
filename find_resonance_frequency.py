# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 2022 at 12:34pm 

@author: Paul Lavender-Jones
"""
# %% Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

%matplotlib widget
#%%
t = 1.12e-3
cE_33 = 60e9
rho = 7500

sigmaE = 0
k = 0.5

f_r = 1/(2*t)*np.sqrt(cE_33/rho*(1-sigmaE)/((1+sigmaE)*(1-2*sigmaE)))
f_a = 1/(2*t)*np.sqrt(cE_33/(rho*(1-k*t**2)))

print(format(f_r,'.2E'))
print(format(f_a,'.2E'))

# %%
c_ = cE_33
n = 2
R = 0
k_iq = 0
f_1 = 1/(2*t)*(1-R-(4*k_iq**2)/np.pi**2)*np.sqrt(c_**n/rho)

print(format(f_1,'.2E'))

# %%
