#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 15:25:42 2021

@author: david
"""

#%% Load packages
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% Calculation of strat-trop flux from Boering2004

# VSMOW isotope ratios
R18SMOW = 0.0020052
R17SMOW = 0.0003799

# Mixing ratios of CO2 and N2O in ppmv
mrCO2 = 380
mrN2O = 0.3

# Lambda values from Boering2004 and Young2014
lmB = 0.516
lmY = 0.528

# Isotopic ratio rquations

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
    return D17

# ER-2 measurements from Boering2004

# Altitude in kilometers
hB = np.array([19.1, 19.2, 19.2, 19.2, 19.2, 19.1, 15.6, 17.7, 19.6, 16.8,
               19.4, 19.4, 13.7, 11.3, 18.6, 16.4, 18.9, 18.4, 19.7, 20.5,
               20.8, 14.5, 18.6, 20.6, 19.3, 19.4, 19.2, 15.8, 20])

# N2O mixing ratio in ppb
N2OB = np.array([154, 187, 249, 205, 128, 97, 217, 141, 131, 263, 198, 218,
                292, 303, 195, 237, 170, 223, 220, 239, 104, 296, 260, 82,
                211, 205, 195, 274, 210])

# d'17O of stratospheric CO2 relative to VSMOW
d17CO2sB = np.array([25.13, 23.70, 22.28, 22.93, 24.31, 26.19, 23.46, 26.28,
                    25.31, 22.28, 23.92, 23.71, 20.58, 21.11, 22.86, 21.74,
                    23.78, 23.15, 23.05, 23.06, 28.20, 19.24, 22.31, 28.61,
                    23.66, 23.37, 24.62, 21.83, 22.93])
dp17CO2sB = dXdpX(d17CO2sB)

# d'18O of stratospheric CO2 relative to VSMOW
d18CO2sB = np.array([41.52, 40.27, 40.15, 39.78, 39.14, 40.76, 40.74, 42.71,
                   40.41, 40.63, 41.03, 41.40, 38.57, 39.78, 39.05, 38.54,
                   39.51, 40.84, 40.38, 41.03, 43.75, 36.43, 40.67, 43.57,
                   41.16, 40.34, 41.83, 40.38, 40.05])
dp18CO2sB = dXdpX(d18CO2sB)

# D17 of stratospheric CO2
D17CO2sBB = D17(d17CO2sB, d18CO2sB, lmB)
D17CO2sBY = D17(dp17CO2sB, dp18CO2sB, lmY)

# Compiled Boering data
Bdata = pd.DataFrame({'hB':hB, 'N2OB':N2OB, 'd17CO2sB':d17CO2sB,
                      'dp17CO2sB':dp17CO2sB,'d18CO2sB':d18CO2sB,
                      'dp18CO2sB':dp18CO2sB, 'D17CO2sBB':D17CO2sBB,
                      'D17CO2sBY': D17CO2sBY})
Bdata = Bdata.sort_values('hB', 0).reset_index(drop=True)

# Boering data with N2O > 195
Bdata195 = Bdata[Bdata['N2OB'] > 195]

# Linear regression of stratospheric CO2 D17 vs. N2O mixing ratio
regBB = stats.linregress(Bdata195['N2OB'], Bdata195['D17CO2sBB'])
regBY = stats.linregress(Bdata195['N2OB'], Bdata195['D17CO2sBY'])

# Vertical flux of N2O in units of moles of N
N2Overt = 13e12

# Vertical flux of CO2 in grams of N
CO2vert = N2Overt / (np.absolute(regBB[0]) * (mrN2O / mrCO2))

# Vertical flux of CO2 in moles of N
CO2vert = CO2vert / 44.01

# Vertical flux of D17CO2s

# Boering xJ
xJB = 3.6e15

#%% Calculation of strat-trop flux from Wiegel2013

# Altitude in kilometers
hW = np.array([16.2, 11.2, 19.6, 19.4, 20.1, 19.4, 19.5, 16.2, 20.2, 19.8,
               12.3, 19.2, 19.5, 17.9, 18.3, 19.9, 19.6, 19.7, 19.3, 17.1,
               17.6, 18.7, 14.5, 11.7, 19, 17.1, 16.8, 19.9, 17.2, 19.5,
               19.1, 19.3, 18.5, 18.9, 19.2, 33.3, 32.2, 31.5, 30.8, 30,
               29.2, 28.7, 28, 27.2])

# N2O mixing ratio in ppb
N2OW = np.array([301, 315, 281, 256, 211, 228, 222, 286, 120, 139, 296, 175,
                 149, 195, 180, 142, 145, 129, 150, 208, 193, 171, 293, 305,
                 149, 221, 241, 70, 279, 106, 109, 112, 232, 252, 216, 54,
                 56.7, 79.6, 77.8, 80.2, 92.1, 107.9, 135.2, 155.9])

# d'17O of stratospheric CO2 relative to VSMOW
dp17CO2sW = np.array([21.7, 21.2, 22.7, 23.8, 24.4, 24.2, 24.2, 22, 28, 27.3,
                      22.3, 25.6, 26.8, 24.9, 25.2, 26.4, 26.4, 27.5, 25.7,
                      24.7, 25.5, 25.6, 21.7, 21.4, 26.6, 23.4, 23.7, 29.1,
                      21.9, 27.7, 27.4, 28.2, 23.1, 22.7, 24.2, 29.1, 29.1,
                      28, 28.6, 28.7, 28.7, 27.5, 26.8, 26.5])
d17CO2sW = dpXdX(dp17CO2sW)

# d'18O of stratospheric CO2 relative to VSMOW
dp18CO2sW = np.array([40.59, 40.24, 41.16, 41.64, 41.7, 41.93, 41.98, 40.87,
                      43.61, 41.14, 40.83, 42.24, 43.17, 42.22, 42.31, 42.72,
                      42.81, 43.18, 42.53, 42.07, 42.19, 42.45, 40.86, 40.52,
                      43.02, 41.84, 41.84, 44.27, 41.11, 43.44, 43.26, 43.61,
                      41.87, 41.14, 41.97, 44.17, 44.46, 43.84, 44.18, 44.44,
                      44.02, 43.84, 43.33, 43.17])
d18CO2sW = dpXdX(dp18CO2sW)

# D17 of stratospheric CO2
D17CO2sW = D17(dp17CO2sW, dp18CO2sW, lmY)
D17CO2sWB = D17(d17CO2sW, d18CO2sW, lmB)

# Compiled Wiegel data
Wdata = pd.DataFrame({'hW':hW, 'N2OW':N2OW, 'dp17CO2sW':dp17CO2sW,
                      'd17CO2sW': d17CO2sW, 'dp18CO2sW':dp18CO2sW,
                      'd18CO2sW':d18CO2sW, 'D17CO2sW':D17CO2sW,
                      'D17CO2sWB': D17CO2sWB})
Wdata = Wdata.sort_values('hW', 0).reset_index(drop=True)

# Wiegel data with N2O > 195
Wdata195 = Wdata[Wdata['N2OW'] > 195]

# Linear regression of stratospheric CO2 D17 vs. N2O mixing ratio
regWB = stats.linregress(Wdata195['N2OW'], Wdata195['D17CO2sWB'])
regWY = stats.linregress(Wdata195['N2OW'], Wdata195['D17CO2sW'])

# Wiegel xJ scaled from Boering xJ
xJW = xJB * (regWB[0] / regBB[0])

#%% Calculation of altitude integrated d18CO2t and D17CO2t

## Loading data

# Liang 2017 data
Liang2017 = pd.read_excel('CO2t.xlsx', 0)
Liang2017 = Liang2017[['D17', 'd18O']]

# Barkan Luz 2012 data
BarkanLuz2012 = pd.read_excel('CO2t.xlsx', 1)
BarkanLuz2012 = BarkanLuz2012[['d17O', 'd18O']]

# Thiemens 2014 data
Thiemens2014 = pd.read_excel('CO2t.xlsx', 2)
Thiemens2014 = Thiemens2014[['d17O', 'd18O']]

## Converting d to d'

# Liang 2017 d'18O and d'17O
d17CO2tL = np.exp(Liang2017['D17'] + 0.516 * np.log(Liang2017['d18O'] + 1)) - 1
d17CO2tL = Liang2017['D17'] + Liang2017['d18O'] * 0.516
dp17CO2tL = dXdpX(d17CO2tL)
dp18CO2tL = dXdpX(Liang2017['d18O'])

# Barkan Luz 2012 d'17O and d'18O

# on VSMOW-CO2 scale
d17CO2tBL = BarkanLuz2012['d17O']
d18CO2tBL = BarkanLuz2012['d18O']

# on VSMOW scale
aCO2H2O = 1.0412
d17CO2tBL = (d17CO2tBL + 1) / aCO2H2O - 1
d18CO2tBL = (d18CO2tBL + 1) / aCO2H2O - 1

# d'17O and d'18O
dp17CO2tBL = dXdpX(d17CO2tBL)
dp18CO2tBL = dXdpX(d18CO2tBL)

# Thiemens 2014 d'17O and d'18O
d17CO2tT = Thiemens2014['d17O']
d18CO2tT = Thiemens2014['d18O']
dp17CO2tT = dXdpX(Thiemens2014['d17O'])
dp18CO2tT = dXdpX(Thiemens2014['d18O'])

# Combined data
dp17CO2tA = dp17CO2tL
#dp17CO2tA = dp17CO2tA.append(dp17CO2tBL)
dp17CO2tA = dp17CO2tA.append(dp17CO2tT) 
dp17CO2tAavg = np.average(dp17CO2tA)

dp18CO2tA = dp18CO2tL
#dp18CO2tA = dp18CO2tA.append(dp18CO2tBL)
dp18CO2tA = dp18CO2tA.append(dp18CO2tT) 
dp18CO2tAavg = np.average(dp18CO2tA)

## D'17

# Liang 2017 D'17
Dp17CO2tL = D17(dp17CO2tL, dp18CO2tL, lmY)
Dp17CO2tLavg = np.average(Dp17CO2tL)

# Barkan Luz 2012 D'17
Dp17CO2tBL = D17(dp17CO2tBL, dp18CO2tBL, lmY)
Dp17CO2tBLavg = np.average(Dp17CO2tBL)

# Thiemens 2014 D'17
Dp17CO2tT = D17(dp17CO2tT, dp18CO2tT, lmY)
Dp17CO2tTavg = np.average(Dp17CO2tT)

# Combined D'17
Dp17CO2tA = Dp17CO2tL
#Dp17CO2tA = Dp17CO2tA.append(Dp17CO2tBL)
Dp17CO2tA = Dp17CO2tA.append(Dp17CO2tT) 
Dp17CO2tAavg = np.average(Dp17CO2tA)



#%% Calculation of altitude integrated d18CO2s and D17CO2s

# Constants needed
avo = 6.0221367e23 # Avogadro's number in molecules / mol
volT = 2.8e15 # volume of troposphere in cm^3
pCO2 = 0.0375 # percent of troposphere made up from CO2
mmair = 28.07 # molar mass of air in g / mol
kB = 1.380658e-19 # Boltzmann's constant in cm^3 hPA / K molecule

# NOAA 1976 Standard Atmosphere data
NOAA = pd.read_excel('pressureheight.xlsx')

# Number of molecules from NOAA1976 pressure data
nJ76 = (NOAA['p (hPa)'] / (kB * NOAA['T (K)']) * volT)
NOAA['molecules'] = nJ76

# Boering average d18O and D17 of CO2s

# Splitting Boering data into 8 unequal intervals
binsB1 = Bdata[0:1]
binsB3 = Bdata[1:3]
binsB7 = Bdata[3:7]
binsB11 = Bdata[7:11]
binsB15 = Bdata[11:15]
binsB19 = Bdata[15:19]
binsB24 = Bdata[19:24]
binsB29 = Bdata[24:29]


# Splitting NOAA 1976 data into corresponding altitude intervals
binsNB1 = NOAA[11:12]
binsNB3 = NOAA[14:16]
binsNB7 = NOAA[16:18]
binsNB11 = NOAA[18:19]
binsNB15 = NOAA[19:20]
binsNB19 = NOAA[19:20]
binsNB24 = NOAA[19:21]
binsNB29 = NOAA[20:22]

# Integration of molecules in each bin
intnNB1 = float(binsNB1['molecules'])
intnNB3 = float(np.trapz(binsNB3['molecules'], binsNB3['h (km)']))
intnNB7 = float(np.trapz(binsNB7['molecules'], binsNB7['h (km)']))
intnNB11 = float(binsNB11['molecules'])
intnNB15 = float(binsNB15['molecules'])
intnNB19 = float(binsNB19['molecules'])
intnNB24 = float(np.trapz(binsNB24['molecules'], binsNB24['h (km)']))
intnNB29 = float(np.trapz(binsNB29['molecules'], binsNB29['h (km)']))
intnNOAA = float(intnNB1 + intnNB3 + intnNB7 + intnNB11 + intnNB15 +
                 intnNB19 + intnNB24 + intnNB29)

# Weighting by molecules
wnNB1= intnNB1 / intnNOAA
wnNB3 = intnNB3 / intnNOAA
wnNB7 = intnNB7 / intnNOAA
wnNB11 = intnNB11 / intnNOAA
wnNB15 = intnNB15 / intnNOAA
wnNB19 = intnNB19 / intnNOAA
wnNB24 = intnNB24 / intnNOAA
wnNB29 = intnNB29 / intnNOAA

# Average d18CO2s in each bin
d18CO2savgB1 = float(np.average(binsB1['d18CO2sB']))
d18CO2savgB3 = float(np.average(binsB3['d18CO2sB']))
d18CO2savgB7 = float(np.average(binsB7['d18CO2sB']))
d18CO2savgB11 = float(np.average(binsB11['d18CO2sB']))
d18CO2savgB15 = float(np.average(binsB15['d18CO2sB']))
d18CO2savgB19 = float(np.average(binsB19['d18CO2sB']))
d18CO2savgB24 = float(np.average(binsB24['d18CO2sB']))
d18CO2savgB29 = float(np.average(binsB29['d18CO2sB']))

# Average d'18CO2s in each bin
dp18CO2savgB1 = float(np.average(binsB1['dp18CO2sB']))
dp18CO2savgB3 = float(np.average(binsB3['dp18CO2sB']))
dp18CO2savgB7 = float(np.average(binsB7['dp18CO2sB']))
dp18CO2savgB11 = float(np.average(binsB11['dp18CO2sB']))
dp18CO2savgB15 = float(np.average(binsB15['dp18CO2sB']))
dp18CO2savgB19 = float(np.average(binsB19['dp18CO2sB']))
dp18CO2savgB24 = float(np.average(binsB24['dp18CO2sB']))
dp18CO2savgB29 = float(np.average(binsB29['dp18CO2sB']))

# Average D17CO2s in each bin
D17CO2savgB1 = float(np.average(binsB1['D17CO2sBB']))
D17CO2savgB3 = float(np.average(binsB3['D17CO2sBB']))
D17CO2savgB7 = float(np.average(binsB7['D17CO2sBB']))
D17CO2savgB11 = float(np.average(binsB11['D17CO2sBB']))
D17CO2savgB15 = float(np.average(binsB15['D17CO2sBB']))
D17CO2savgB19 = float(np.average(binsB19['D17CO2sBB']))
D17CO2savgB24 = float(np.average(binsB24['D17CO2sBB']))
D17CO2savgB29 = float(np.average(binsB29['D17CO2sBB']))

# Average D'17CO2s in each bin
Dp17CO2savgB1 = float(np.average(binsB1['D17CO2sBY']))
Dp17CO2savgB3 = float(np.average(binsB3['D17CO2sBY']))
Dp17CO2savgB7 = float(np.average(binsB7['D17CO2sBY']))
Dp17CO2savgB11 = float(np.average(binsB11['D17CO2sBY']))
Dp17CO2savgB15 = float(np.average(binsB15['D17CO2sBY']))
Dp17CO2savgB19 = float(np.average(binsB19['D17CO2sBY']))
Dp17CO2savgB24 = float(np.average(binsB24['D17CO2sBY']))
Dp17CO2savgB29 = float(np.average(binsB29['D17CO2sBY']))

# Weighted average d18CO2s
d18CO2savgB = (d18CO2savgB1 * wnNB1 + d18CO2savgB3 * wnNB3 +
               d18CO2savgB7 * wnNB7 + d18CO2savgB11 * wnNB11 +
               d18CO2savgB15 * wnNB15 + d18CO2savgB19 * wnNB19 + 
               d18CO2savgB24 * wnNB24 + d18CO2savgB29 * wnNB29)
print(d18CO2savgB)

# Weighted average dp18CO2s
dp18CO2savgB = (dp18CO2savgB1 * wnNB1 + dp18CO2savgB3 * wnNB3 +
                dp18CO2savgB7 * wnNB7 + dp18CO2savgB11 * wnNB11 +
                dp18CO2savgB15 * wnNB15 + dp18CO2savgB19 * wnNB19 + 
                dp18CO2savgB24 * wnNB24 + dp18CO2savgB29 * wnNB29)
print(dp18CO2savgB)

# Weighted average D17CO2s
D17CO2savgB = (D17CO2savgB1 * wnNB1 + D17CO2savgB3 * wnNB3 +
               D17CO2savgB7 * wnNB7 + D17CO2savgB11 * wnNB11 +
               D17CO2savgB15 * wnNB15 + D17CO2savgB19 * wnNB19 +
               D17CO2savgB24 * wnNB24 + D17CO2savgB29 * wnNB29)
print(D17CO2savgB)

# Weighted average Dp17CO2s
Dp17CO2savgB = (Dp17CO2savgB1 * wnNB1 + Dp17CO2savgB3 * wnNB3 +
                Dp17CO2savgB7 * wnNB7 + Dp17CO2savgB11 * wnNB11 +
                Dp17CO2savgB15 * wnNB15 + Dp17CO2savgB19 * wnNB19 +
                Dp17CO2savgB24 * wnNB24 + Dp17CO2savgB29 * wnNB29)
print(Dp17CO2savgB)

# Weighted average d18CO2s
Dp17CO2savgB1 = stats.binned_statistic(hB, D17CO2sBY)

# Model output in Boering reference frame
dp18CO2s = 43.78
d18CO2s = dpXdX(dp18CO2s)
dp17CO2s = 25.49
d17CO2s = dpXdX(dp17CO2s)
Dp17CO2s = D17(dp17CO2s, dp18CO2s, lmY)
D17CO2s = D17(d17CO2s, d18CO2s, lmB)

# Wiegel average d18O and D17 of CO2s

# Splitting Wiegel data into 8 unequal intervals
binsW6 = Wdata[0:6]
binsW12 = Wdata[6:12]
binsW18 = Wdata[12:18]
binsW24 = Wdata[18:24]
binsW30 = Wdata[24:30]
binsW35 = Wdata[30:35]
binsW40 = Wdata[35:40]
binsW44 = Wdata[40:44]

# Splitting NOAA1976 data into corresponding altitude intervals
binsNW6 = NOAA[11:17]
binsNW12 = NOAA[17:19]
binsNW18 = NOAA[18:20]
binsNW24 = NOAA[19:20]
binsNW30 = NOAA[19:21]
binsNW35 = NOAA[20:21]
binsNW40 = NOAA[27:31]
binsNW44 = NOAA[31:34]

# Integration of molecules in each bin
intnNW6 = float(np.trapz(binsNW6['molecules'], binsNW6['h (km)']))
intnNW12 = float(np.trapz(binsNW12['molecules'], binsNW12['h (km)']))
intnNW18 = float(np.trapz(binsNW18['molecules'], binsNW18['h (km)']))
intnNW24 = float(binsNW24['molecules'])
intnNW30 = float(np.trapz(binsNW30['molecules'], binsNW30['h (km)']))
intnNW35 = float(binsNW35['molecules'])
intnNW40 = float(np.trapz(binsNW40['molecules'], binsNW40['h (km)']))
intnNW44 = float(np.trapz(binsNW44['molecules'], binsNW44['h (km)']))
intnNOAA = float(intnNW6 + intnNW12 + intnNW18 + intnNW24 + intnNW30 +
                 intnNW35 + intnNW40 + intnNW44)

# Weighting by molecules
wnNW6 = intnNW6 / intnNOAA
wnNW12 = intnNW12 / intnNOAA
wnNW18 = intnNW18 / intnNOAA
wnNW24 = intnNW24 / intnNOAA
wnNW30 = intnNW30 / intnNOAA
wnNW35 = intnNW35 / intnNOAA
wnNW40 = intnNW40 / intnNOAA
wnNW44 = intnNW44 / intnNOAA

# Average d18CO2s in each bin
d18CO2savgW6 = float(np.average(binsW6['d18CO2sW']))
d18CO2savgW12 = float(np.average(binsW12['d18CO2sW']))
d18CO2savgW18 = float(np.average(binsW18['d18CO2sW']))
d18CO2savgW24 = float(np.average(binsW24['d18CO2sW']))
d18CO2savgW30 = float(np.average(binsW30['d18CO2sW']))
d18CO2savgW35 = float(np.average(binsW35['d18CO2sW']))
d18CO2savgW40 = float(np.average(binsW40['d18CO2sW']))
d18CO2savgW44 = float(np.average(binsW44['d18CO2sW']))

# Average d'18CO2s in each bin
dp18CO2savgW6 = float(np.average(binsW6['dp18CO2sW']))
dp18CO2savgW12 = float(np.average(binsW12['dp18CO2sW']))
dp18CO2savgW18 = float(np.average(binsW18['dp18CO2sW']))
dp18CO2savgW24 = float(np.average(binsW24['dp18CO2sW']))
dp18CO2savgW30 = float(np.average(binsW30['dp18CO2sW']))
dp18CO2savgW35 = float(np.average(binsW35['dp18CO2sW']))
dp18CO2savgW40 = float(np.average(binsW40['dp18CO2sW']))
dp18CO2savgW44 = float(np.average(binsW44['dp18CO2sW']))

# Average D17CO2s in each bin
D17CO2savgW6 = float(np.average(binsW6['D17CO2sWB']))
D17CO2savgW12 = float(np.average(binsW12['D17CO2sWB']))
D17CO2savgW18 = float(np.average(binsW18['D17CO2sWB']))
D17CO2savgW24 = float(np.average(binsW24['D17CO2sWB']))
D17CO2savgW30 = float(np.average(binsW30['D17CO2sWB']))
D17CO2savgW35 = float(np.average(binsW35['D17CO2sWB']))
D17CO2savgW40 = float(np.average(binsW40['D17CO2sWB']))
D17CO2savgW44 = float(np.average(binsW44['D17CO2sWB']))

# Average D'17CO2s in each bin
Dp17CO2savgW6 = float(np.average(binsW6['D17CO2sW']))
Dp17CO2savgW12 = float(np.average(binsW12['D17CO2sW']))
Dp17CO2savgW18 = float(np.average(binsW18['D17CO2sW']))
Dp17CO2savgW24 = float(np.average(binsW24['D17CO2sW']))
Dp17CO2savgW30 = float(np.average(binsW30['D17CO2sW']))
Dp17CO2savgW35 = float(np.average(binsW35['D17CO2sW']))
Dp17CO2savgW40 = float(np.average(binsW40['D17CO2sW']))
Dp17CO2savgW44 = float(np.average(binsW44['D17CO2sW']))

# Weighted average d18CO2s
d18CO2savgW = (d18CO2savgW6 * wnNW6 + d18CO2savgW12 * wnNW12 +
               d18CO2savgW18 * wnNW18 + d18CO2savgW24 * wnNW24 +
               d18CO2savgW30 * wnNW30 + d18CO2savgW35 * wnNW35 +
               d18CO2savgW40 * wnNW40 + d18CO2savgW44 * wnNW44)
print(d18CO2savgW)

# Weighted average d18CO2s
dp18CO2savgW = (dp18CO2savgW6 * wnNW6 + dp18CO2savgW12 * wnNW12 +
                dp18CO2savgW18 * wnNW18 + dp18CO2savgW24 * wnNW24 +
                dp18CO2savgW30 * wnNW30 + dp18CO2savgW35 * wnNW35 +
                dp18CO2savgW40 * wnNW40 + dp18CO2savgW44 * wnNW44)
print(dp18CO2savgW)

# Weighted average D17CO2s
D17CO2savgW = (D17CO2savgW6 * wnNW6 + D17CO2savgW12 * wnNW12 +
               D17CO2savgW18 * wnNW18 + D17CO2savgW24 * wnNW24 +
               D17CO2savgW30 * wnNW30 + D17CO2savgW35 * wnNW35 +
               D17CO2savgW40 * wnNW40 + D17CO2savgW44 * wnNW44)
print(D17CO2savgW)

# Weighted average Dp17CO2s
Dp17CO2savgW = (Dp17CO2savgW6 * wnNW6 + Dp17CO2savgW12 * wnNW12 +
                Dp17CO2savgW18 * wnNW24 + Dp17CO2savgW24 * wnNW24 +
                Dp17CO2savgW30 * wnNW30 + Dp17CO2savgW35 * wnNW35 +
                Dp17CO2savgW40 * wnNW40 + Dp17CO2savgW44 * wnNW44 )
print(Dp17CO2savgW)