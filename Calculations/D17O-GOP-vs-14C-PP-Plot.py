#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 09:51:31 2021

This file plots estimates of gross oxygen production based on oxygen triple-isotope measurements of seawater against estimates
of net carbon fixation based on C14 measurements. It is Figure 10 in Valerio 2021.

@author: david
"""

#%% Load packages
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines

#%% D17O GOP and 14C-PP data to plot 

# Luz and Barkan 2009 BATS
BATSGOP = np.average(np.array([37.1, 47.5, 76.2]))
BATSC14 = np.average(np.array([8.6, 7.6, 9.9]))

# Munro et al. 2012 CalCOFI
CalCOFIGOP = 235
CalCOFIC14 = 529 * (1 / 1000) * (1/12) * 1000

# Stanley et al. 2010 Central Equatorial Pacific (CEP)
CEPGOP = 184
CEPGOPC14 = 3.1
CEPC14 = CEPGOP / CEPGOPC14

# Stanley et al. 2010 Western Equatorial Pacific (WEP)
WEPGOP = 121
WEPGOPC14 = 8.2
WEPC14 = WEPGOP / WEPGOPC14

# Quay 2010 HOT
HOTGOP = 103
HOTC14 = 42

# Hamme et al. 2012 Southern Ocean (Patch 1 + 2)
SOGOP = np.average(np.array([144, 159]))
SOC14 = np.average(np.array([39.1, 27]))

#%% GOP / C14 slopes to plot

# JGOFS and Halsey2010 GOPC14 slopes
Halsey2010 = 3.3
JGOFS = 2.7

# Slopes with 50% Mehler consumption
Halsey2010MR = Halsey2010 * 1.5
JGOFSMR = JGOFS * 1.5

#%% Plots

# Setting up figure parameters
fig1= plt.figure(figsize = (5, 5))
fig1 = fig1.add_subplot(1, 1, 1)
fig1.set(xlim = (0, 100), ylim = (0, 500))
fig1.set_facecolor('#F7F7F7')

#fig1.grid()

# Plot GOP / C14 slopes
x = np.linspace(0, 100)

fig1.plot(x, Halsey2010 * x, color='#a50026', label='Halsey2010', zorder=3.5)
fig1.text(63, 223, 'Halsey2010', color='#a50026', fontsize=8,
          rotation=36)
fig1.plot(x, Halsey2010MR * x, color='#a50026', linestyle='dashed',
          label='Halsey2010 + 50%', zorder=3.5)
fig1.text(57, 305, 'Halsey2010 + 50%', color='#a50026', fontsize=8,
          rotation=46)
fig1.plot(x, JGOFS * x, color='#364b9a', label='JGOFS', zorder=3.5)
fig1.text(65,  189, 'JGOFS', color='#364b9a', fontsize=8,
          rotation=30)
fig1.plot(x, JGOFSMR * x, color='#364b9a', linestyle='dashed',
          label='JGOFS + 50%', zorder=3.5)
fig1.text(60,  262, 'JGOFS + 50%', color='#364b9a', fontsize=8,
          rotation=40)

# Plot D17O GOP and 14C-PP data
fig1.scatter(BATSC14, BATSGOP, color='black', label='BATS', zorder=4.5)
fig1.text(BATSC14 - 7, BATSGOP + 15, 'BATS', fontsize=8,
          color='black', zorder=4.5)
fig1.scatter(CalCOFIC14, CalCOFIGOP, color='black', marker='^', label='CalCOFI',
             zorder=4.5)
fig1.text(CalCOFIC14 - 13, CalCOFIGOP + 5, 'CalCOFI', fontsize=8,
          color='black', zorder=4.5)
fig1.scatter(CEPC14, CEPGOP, color='black', marker='s', label='CEP',
             zorder=4.5)
fig1.text(CEPC14 + 5, CEPGOP - 30, 'Central EP', fontsize=8,
          color='black', zorder=4.5)
fig1.scatter(WEPC14, WEPGOP, color='black', marker='p', label='WEP',
             zorder=4.5)
fig1.text(WEPC14 - 12, WEPGOP + 15, 'Western EP', fontsize=8,
          color='black', zorder=4.5)
fig1.scatter(HOTC14, HOTGOP, color='black', marker='h', label='HOTS',
             zorder=4.5)
fig1.text(HOTC14 + 2, HOTGOP - 20, 'HOT', fontsize=8,
          color='black', zorder=4.5)
fig1.scatter(SOC14, SOGOP, color='black', marker='D', label='SO',
             zorder=4.5)
fig1.text(SOC14 - 25, SOGOP + 20, 'Southern Ocean', fontsize=8,
          color='black', zorder=4.5)

# Title and x and y labels
fig1.set_xlabel('Mixed-layer integrated $^{14}$C-PP (mmol C m$^{-2}$ day$^{-1}$)')
fig1.set_ylabel('Mixed-layer $^{17}\Delta$-GOP (mmol O$_{2}$ m$^{-2}$ day$^{-2}$)')

# Legend

# Labels
#labelLines(plt.gca().get_lines(), zorder=3.5, fontsize=8,
           #xvals=[78, 73, 80, 78])
plt.tight_layout()

# Saving figure
plt.savefig('O2C.jpg', dpi=800)



