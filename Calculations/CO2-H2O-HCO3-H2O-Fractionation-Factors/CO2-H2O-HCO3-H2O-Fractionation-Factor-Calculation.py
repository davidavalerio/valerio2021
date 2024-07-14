#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:58:04 2021

This file defines functions for calculating the fractionaction factor between CO2 and H2O and HCO3- and H20 as a function
of temperature and plots them as a function of temperature.

@author: david
"""

#%% Load packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#%% Defining functions for equations calculating CO2-H2O and HCO3-H2O alphas as a function of temperature from Beck et al. 2005

def aCO2H2O(T):
    aCO2H2O = np.exp((2.52 * (10 ** 6 * T ** -2) + 12.12) / 1000)
    return aCO2H2O;
    
def aHCO3H2O(T):
    aHCO2H2O = np.exp((2.59 * (10 ** 6 * T ** -2) + 1.89) / 1000)
    return aCO2H2O

#%% Plotting the fractionation factors as a function of temperature

TC = np.linspace(5, 40)
TK = TC + 273.15

aCO2H2O = aCO2H2O(TK)
aHCO3H2O = aHCO3H2O(TK)

# aCO2-H2O as a function of temperature

# Setting up figure parameters
fig1 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig1 = fig1.add_subplot(1, 1, 1)
fig1.set(xlim = (5, 40), ylim = (1.035, 1.05))

# Plotting
fig1.plot(TC, aCO2H2O)

# Legend and title
fig1.set_xlabel("Temperature (°C)")
fig1.set_ylabel("α$_{CO2-H2O}$")
fig1.legend(loc = 'best', fontsize = 'small')

# Saving plot
plt.savefig('aCO2H2O.png', dpi = 800)

# aHCO3-H2O as a function of temperature

# Setting up figure parameters
fig1 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig1 = fig1.add_subplot(1, 1, 1)
fig1.set(xlim = (5, 40), ylim = (1.035, 1.05))

# Plotting
fig1.plot(TC, aCO2H2O)

# Legend and title
fig1.set_xlabel("Temperature (°C)")
fig1.set_ylabel("α$_{CO2-H2O}$")
fig1.legend(loc = 'best', fontsize = 'small')

# Saving plot
plt.savefig('aCO2H2O.png', dpi = 800)
