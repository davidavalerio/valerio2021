#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 10:27:33 2020

This file produces a contour plot of the oxygen triple-isotope compositions of tropospheric oxygen for different combinations of
global average respiration and evapotranspiration theta values.

@author: david
"""

# Load packages
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# D17 O2 troposphere as function of thetas linear regressin parameters
D17tCO2H2Om = 0.1751257079694242
D17tCO2H2Ob = -0.5173411443774644

D17tevapm = 0.02266821781617719
D17tevapb = -0.436662043606977

D17trespm = 23.16944176869028
D17trespb = -12.354638844550406

# Thetas of interest
tCO2H2O = np.linspace(0.505, 0.530)
tevap = np.linspace(0.505, 0.530)
tresp = np.linspace(0.505, 0.530)

# Contour plot equation
def z(t1, t2, m1, m2, b1, b2):
    z = (m1 * t1 - m2 * t2) + (b1 - b2)
    return z

# Contour plot of tCO2H2O vs. tevap

# theta values to evaluate in contour plot
tCO2H2O, tevap = np.meshgrid(tCO2H2O, tevap)

# Setting up figure parameters
fig1 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig1 = fig1.add_subplot(1, 1, 1)
fig1.set(xlim = (.505, .530), ylim = (0.505, 0.530))
fig1.contour(tCO2H2O, tevap, z(tCO2H2O, tevap, D17tCO2H2Om, D17tevapm,
                              D17tCO2H2Ob, D17tevapb))
fig1.clabel(fig1.contour(tCO2H2O, tevap, z(tCO2H2O, tevap, D17tCO2H2Om, D17tevapm,
                              D17tCO2H2Ob, D17tevapb)), inline=False)

# Legend and title
fig1.set_xlabel("$\Theta_{CO2-H2O}$")
fig1.set_ylabel("$\Theta_{evap}$")
plt.tight_layout()

# Export plot
plt.savefig('tCO2H2Ovtevap.jpg', dpi = 800)

# Contour plot of tresp vs. tCO2H2O

# Theta values to evaluate in contour plot
tresp, tCO2H2O = np.meshgrid(tresp, tCO2H2O)

# Setting up figure parameters
fig2 = plt.figure(figsize=(5, 5))
with sns.axes_style("whitegrid"):
    fig2 = fig2.add_subplot(1, 1, 1)
fig2.set(xlim = (.505, .530), ylim = (0.505, 0.530))
fig2.contour(tresp, tCO2H2O, z(tresp, tCO2H2O, D17trespm, D17tCO2H2Om,
                              D17trespb, D17tCO2H2Ob))
fig2.clabel(fig2.contour(tresp, tCO2H2O, z(tresp, tCO2H2O, D17trespm, D17tCO2H2Om,
                              D17trespb, D17tCO2H2Ob)), inline=False)

# Legend and title
fig2.set_xlabel("$\Theta_{resp}$")
fig2.set_ylabel("$\Theta_{CO2-H2O}$")
plt.tight_layout()

# Export plot
plt.savefig('trespvtCO2H2O.jpg', dpi = 800)

# Contour plot of tresp vs. tevap

# Theta values to evaluate in contour plot
tresp, tevap = np.meshgrid(tresp, tevap)

# Setting up figure parameters
fig3 = plt.figure(figsize=(5, 5))
with sns.axes_style("whitegrid"):
    fig3 = fig3.add_subplot(1, 1, 1)
fig3.set(xlim = (.505, .530), ylim = (0.505, 0.530))
fig3.contour(tresp, tevap, z(tresp, tevap, D17trespm, D17tevapm,
                              D17trespb, D17tevapb))
fig3.clabel(fig3.contour(tresp, tevap, z(tresp, tevap, D17trespm, D17tevapm,
                                         D17trespb, D17tevapb)), inline=False)

# Legend and title
fig3.set_xlabel("$\Theta_{resp}$")
fig3.set_ylabel("$\Theta_{evap}$")
plt.tight_layout()

#Export plots
plt.savefig('trespvtevap.jpg', dpi = 800)

