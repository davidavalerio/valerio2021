#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 13:10:46 2021

@author: david
"""

#%% Load packages
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% Data

# NOAA 1976 altitude, temperature, and pressure data
NOAA1976 = pd.read_excel('NOAA1976.xlsx')

# Boering 2004 altitude, N2O, d17O, and d18O data
Boering2004 = pd.read_excel('CO2s.xlsx', 0)

# Wiegel 2013 altitude, N2O, d17O, and d18O data
Wiegel2013 = pd.read_excel('CO2s.xlsx', 1)

# Liang2017 CO2, d18O, and D17 data
Liang2017 = pd.read_excel('CO2t.xlsx', 0)

# BarkanLuz2012 CO2, d17O, d18O, and D17 data
BarkanLuz2012 = pd.read_excel('CO2t.xlsx', 1)

# Thiemens 2013 
Thiemens2013 = pd.read_excel('CO2t.xlsx', 2)

#%% Equations for calculation

# dXO to d'XO
def dXdpX(dX):
    dpX = 1e3 * np.log(dX / 1e3 + 1)
    return dpX;

# d'XO to dXO
def dpXdX(dpX):
    dpX = 1e3 * (np.exp(dpX / 1e3) - 1)
    return dpX;

# D17
def D17(d17, d18, lm):
    D17 = d17 - lm * d18
    return D17;

def d17OD17(D17, lm, d18O):
    d17O = D17 + lm * d18O
    return d17O;

# dX VSMOw-CO2 to dX VSMOW
def dXref(dXo, a):
    dX = a * dXo + (a - 1) * 1000
    return dX;

# Number density from pressure and temperature data
def nd(p, T):
    nd = p / (kB * T)
    return nd;

# Constants needed
lm = 0.528
avo = 6.0221367e23 # Avogadro's number in molecules / mol
volT = 2.8e15 # volume of troposphere in cm^3
pCO2 = 0.0375 # percent of troposphere made up from CO2
mmair = 28.07 # molar mass of air in g / mol
kB = 1.380658e-19 # Boltzmann's constant in cm^3 hPA / K molecule
aCO2H2O = 1.0412 # alpha-18 between CO2 and H2O

#%% d18O, d'18O, d17O, d'17O data

## Stratosphere

# Boering 2004 d18O
d18OB = Boering2004['d18O']

# Boering 2004 d'18O
dp18OB = dXdpX(d18OB)

# Boering 2004 d17O
d17OB = Boering2004['d17O']

# Boering 2004 d'17O
dp17OB = dXdpX(d17OB)

# Wiegel 2013 d18O
d18OW = Wiegel2013['d18O']

# Wiegel 2013 d'18O
dp18OW = dXdpX(d18OW)

# Wiegel 2013 d17O
d17OW = Wiegel2013['d17O']

# Wiegel 2013 d'17O
dp17OW = dXdpX(d17OW)

## Troposphere

# Liang 2017 d18O
d18OL = Liang2017['d18O']

# Liang 2017 d'18O
dp18OL = dXdpX(d18OL)

# Liang 2017 d17O
D17L = Liang2017['D17']
d17OL = d17OD17(D17L, 0.528, d18OL)

# Liang 2017 d'17O
dp17OL = dXdpX(d17OL)

# Barkan and Luz 2012 d18O relative to VSMOW-CO2
d18OBL = BarkanLuz2012['d18O']

# Barkan and Luz 2012 d18O relative to VSMOW
d18OBL = dXref(d18OBL, aCO2H2O)

# Barkan and Luz 2012 d'18O
dp18OBL = dXdpX(d18OBL)

# Barkan and Luz 2012 d17O relative to VSMOW-CO2
d17OBL = BarkanLuz2012['d17O']

# Barkan and Luz 2012 d17O relative to VSMOW
d17OBL = dXref(d17OBL, aCO2H2O ** 0.5)

# Barkan and Luz 2012 d'17O
dp17OBL = dXdpX(d17OBL)

# Thiemens 2013 d18O
d18OT = Thiemens2013['d18O']

# Thiemens 2013 d'18O
dp18OT = dXdpX(d18OT)

# Thiemens 2013 d17O
d17OT = Thiemens2013['d17O']

# Thiemens 2013 d'17O
dp17OT = dXdpX(d17OT)

#%% D17 in correct reference frame

## Stratosphere

# Boering 2004 D'17
Dp17B = D17(dp17OB, dp18OB, lm)

# Wiegel 2013 D'17
Dp17W = D17(dp17OW, dp18OW, lm)

## Troposphere

# Liang 2017 D'17
Dp17L = D17(dp17OL, dp18OL, lm)

# Barkan and Luz 2012 D'17
Dp17BL = D17(dp17OBL, dp18OBL, lm)

# Thiemens 2013 D'17
Dp17T = D17(dp17OT, dp18OT, lm)

#%% Average d'17O, d'18O, and D'17 of CO2 in the troposphere

# Combined d'17O
dp17Otall = dp17OL.append(dp17OBL)
dp17Otall = dp17Otall.append(dp17OT)

# Average d'17O
dp17Otavg = np.average(dp17Otall)

# Combined d'18O
dp18Otall = dp18OL.append(dp18OBL)
dp18Otall = dp18Otall.append(dp18OT)

# Average d'18O
dp18Otavg = np.average(dp18Otall)

# Combined D'17
Dp17tall = Dp17L.append(Dp17BL)
Dp17tall = Dp17tall.append(Dp17T)

# Average D'17
Dp17tavg = np.average(Dp17tall)

#%% Weighted average d'17O, d'18O, and D'17 of CO2 in the stratosphere

## Boering 2004

# Height
hB = Boering2004['height']

# Bin averages for d'18O
dp18OBavg = stats.binned_statistic(hB, dp18OB, 'mean', bins=8)

# Bin averages for d'17O
dp17OBavg = stats.binned_statistic(hB, dp17OB, 'mean', bins=8)

# Bin averages for D'17
Dp17Bavg = stats.binned_statistic(hB, Dp17B, 'mean', bins=8)

## Wiegel 2013

# Height
hW = Wiegel2013['height']

# Bin averages for d'18O
dp18OWavg = stats.binned_statistic(hW, dp18OW, 'mean', bins=8)

# Bin averages for d'17O
dp17OWavg = stats.binned_statistic(hW, dp17OW, 'mean', bins=8)

# Bin averages for D'17
Dp17Wavg = stats. binned_statistic(hW, Dp17W, 'mean', bins=8)

#%% Weighting by number density

## Number density data from NOAA 1976

# Altitude, pressure, and temperature data
hNOAA = NOAA1976['h (km)']
pNOAA = NOAA1976['p (hPa)']
TNOAA = NOAA1976['T (K)']

# Number density as a function of altitude
ndNOAA = nd(pNOAA, TNOAA)

# Number density data with altitude as index
ndNOAA = pd.DataFrame({'h': hNOAA,'nd': ndNOAA}).set_index('h')

## Boering 2004

# Bin edges
binsB = Dp17Bavg[1].round(1)

# Number density in corresponding bins
binB1 = ndNOAA[binsB[0]:binsB[1]]
#binB2 = ndNOAA[binsB[1]:binsB[2]]
binB3 = ndNOAA[binsB[2]:binsB[3]]
binB4 = ndNOAA[binsB[3]:binsB[4]]
binB5 = ndNOAA[binsB[4]:binsB[5]]
binB6 = ndNOAA[binsB[5]:binsB[6]]
binB7 = ndNOAA[binsB[6]:binsB[7]]
binB8 = ndNOAA[binsB[7]:binsB[8]]

# Integrated number density in each bin
intB1 = np.trapz(binB1['nd'], binB1.index)
#intB2 = np.trapz(binB2['nd'], binB2.index)
intB3 = np.trapz(binB3['nd'], binB3.index)
intB4 = np.trapz(binB4['nd'], binB4.index)
intB5 = np.trapz(binB5['nd'], binB5.index)
intB6 = np.trapz(binB6['nd'], binB6.index)
intB7 = np.trapz(binB7['nd'], binB7.index)
intB8 = np.trapz(binB8['nd'], binB8.index)
intB = (intB1 + intB3 + intB4 + intB5 + intB6 + intB7 + intB8)

# Weighted integrated number density in each bin
intB1w = intB1 / intB
#intB2w = intB2 / intB
intB3w = intB3 / intB
intB4w = intB4 / intB
intB5w = intB5 / intB
intB6w = intB6 / intB
intB7w = intB7 / intB
intB8w = intB8 / intB

# Weighted average d'18O
dp18OBwavg = (intB1w * dp18OBavg[0][0] +
              intB3w * dp18OBavg[0][2] + intB4w * dp18OBavg[0][3] +
              intB5w * dp18OBavg[0][4] + intB6w * dp18OBavg[0][5] +
              intB7w * dp18OBavg[0][6] + intB8w * dp18OBavg[0][7])

# Weighted average d'17O
dp17OBwavg = (intB1w * dp17OBavg[0][0] +
              intB3w * dp17OBavg[0][2] + intB4w * dp17OBavg[0][3] +
              intB5w * dp17OBavg[0][4] + intB6w * dp17OBavg[0][5] +
              intB7w * dp17OBavg[0][6] + intB8w * dp17OBavg[0][7])

# Weighted average D'17
Dp17Bwavg = (intB1w * Dp17Bavg[0][0] +
             intB3w * Dp17Bavg[0][2] + intB4w * Dp17Bavg[0][3] +
             intB5w * Dp17Bavg[0][4] + intB6w * Dp17Bavg[0][5] +
             intB7w * Dp17Bavg[0][6] + intB8w * Dp17Bavg[0][7])

## Wiegel 2013

# Bin edges
binsW = Dp17Wavg[1].round(1)

# Number density in corresponding bins
binW1 = ndNOAA[binsW[0]:binsW[1]]
binW2 = ndNOAA[binsW[1]:binsW[2]]
binW3 = ndNOAA[binsW[2]:binsW[3]]
binW4 = ndNOAA[binsW[3]:binsW[4]]
#binW5 = ndNOAA[binsW[4]:binsW[5]]
binW6 = ndNOAA[binsW[5]:binsW[6]]
binW7 = ndNOAA[binsW[6]:binsW[7]]
binW8 = ndNOAA[binsW[7]:binsW[8]]

# Integrated number density in each bin
intW1 = np.trapz(binW1['nd'], binW1.index)
intW2 = np.trapz(binW2['nd'], binW2.index)
intW3 = np.trapz(binW3['nd'], binW3.index)
intW4 = np.trapz(binW4['nd'], binW4.index)
#intW5 = np.trapz(binW5['nd'], binW5.index)
intW6 = np.trapz(binW6['nd'], binW6.index)
intW7 = np.trapz(binW7['nd'], binW7.index)
intW8 = np.trapz(binW8['nd'], binW8.index)
intW = (intW1 + intW2 + intW3 + intW4 + intW6 + intW7 + intW8)

# Weighted integrated number density in each bin
intW1w = intW1 / intW
intW2w = intW2 / intW
intW3w = intW3 / intW
intW4w = intW4 / intW
#intW5w = intW5 / intW
intW6w = intW6 / intW
intW7w = intW7 / intW
intW8w = intW8 / intW

# Weighted average d'18O
dp18OWwavg = (intW1w * dp18OWavg[0][0] + intW2w * dp18OWavg[0][1] +
              intW3w * dp18OWavg[0][2] + intW4w * dp18OWavg[0][3] +
                                         intW6w * dp18OWavg[0][5] +
              intW7w * dp18OWavg[0][6] + intW8w * dp18OWavg[0][7])

# Weighted average d'17O
dp17OWwavg = (intW1w * dp17OWavg[0][0] + intW2w * dp17OWavg[0][1] +
              intW3w * dp17OWavg[0][2] + intW4w * dp17OWavg[0][3] +
                                         intW6w * dp17OWavg[0][5] +
              intW7w * dp17OWavg[0][6] + intW8w * dp17OWavg[0][7])

# Weighted average D'17
Dp17Wwavg = (intW1w * Dp17Wavg[0][0] + intW2w * Dp17Wavg[0][1] +
             intW3w * Dp17Wavg[0][2] + intW4w * Dp17Wavg[0][3] +
                                       intW6w * Dp17Wavg[0][5] +
             intW7w * Dp17Wavg[0][6] + intW8w * Dp17Wavg[0][7])

#%% Collated estimates
isocompCO2 = pd.DataFrame([dp18Otavg, dp17Otavg, Dp17tavg, dp18OWwavg,
                           dp17OWwavg, Dp17Wwavg],
                          index=['d18_CO2t', 'd17_CO2t', 'D17_CO2t',
                                 'd18_CO2s', 'd17_CO2s', 'D17_CO2s'])
print(isocompCO2)

#%% xJ calculation

# Wiegel 2013 xJ using regular delta and slope = 0.528 (less than `% diff)
xJW = 49e15 # per mil grams C
xJW = xJW / 12.011 # per mil moles C

