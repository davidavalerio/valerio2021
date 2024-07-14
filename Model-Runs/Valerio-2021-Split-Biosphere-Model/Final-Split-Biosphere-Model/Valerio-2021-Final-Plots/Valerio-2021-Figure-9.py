#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 12:11:16 2021

This file runs Valerio 2021's split biosphere box model to plot how the oxygen-triple isotope composition of the troposphere changes
as a function of different contributions of Mehler-like reactions to global oxygen uptake. It is Figure 9 in Valerio 2021.

@author: david
"""

#%% Load packages
import pandas as pd
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#%% Terrestrial biosphere

# Carbon flux from terrestrial biosphere
FCt = 9.4e15 # mol C / yr from Zhang1996

# Relative fraction of terrestrial productivity for C3 and C4 plants
fC3 = 0.725 # fraction from Francois1998
fC4 = 0.275 # fraction from Francois1998

# Photosynthetic quotient
pQ = 1.07 # mol O2 / mol C from Keeling 1988

# Percentage of O2 production taken by Mehler reaction in terrestrial biosphere
fMRt = 0.1 # from Badger2000 

# AOX relative to COX
fAOXCOX = 0.1 # from Angert2003

# C4 O2 production
FO2C4t = FCt * pQ * fC4

# O2 production in C4 plants taken by Mehler reaction
FO2C4MRt = FO2C4t * (1 / (1 - fMRt)) - FO2C4t

# (dark respiration + photorespiration) / (dark respiration)
def DRPR(CO2, pipa, T):
    DRPR = (4.5 / 4) * (CO2 * pipa + (7 / 3) * T) / (CO2 * pipa - T)
    return DRPR;

# Parameters for DRPR
CO2 = 281 # pre-industrial ppm CO2 in the atmosphere
pipa = 0.65 # ratio of partial pressure CO2 in leaf vs. atmosphere
T = 34 # CO2 compensation point at 28 oC

DRPR = DRPR(CO2, pipa, T)

# O2 production in C3 plants taken by photorespiration
FO2C3t = FCt * pQ * fC3 * DRPR

# O2 production in C3 plants taken by Mehler reaction
FO2C3MRt = FO2C3t * (1 / (1 - fMRt)) - FO2C3t

# O2 production in terrestrial biosphere taken by Mehler reaction
FO2MRt = FO2C4MRt + FO2C3MRt

# O2 production in terrestrial biosphere taken by photorespiration
FO2PRt = FO2C3t - (FO2C3t / DRPR)

# Total terrestial biosphere O2 production
FO2t = FO2C4t + FO2C3t + FO2MRt

#%% Marine biosphere

# Carbon flux from marine biosphere
FCm = 4.04e15 # mol C / yr from Field1998

# Ratio of 14C production to gross O2 production
C14O2 = 0.37 # from Blunier 2002

# Percentage of O2 production taken by Mehler reaciton in marine biosphere
fMRm = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5]) # from Badger2000

# Marine biosphere production w/o MR
FO2m = FCm * (1 / C14O2)

# O2 production in marine biosphere taken by Mehler reaction
FO2MRm0 = FO2m * (1 / (1 - fMRm[0])) - FO2m 
FO2MRm1 = FO2m * (1 / (1 - fMRm[1])) - FO2m 
FO2MRm2 = FO2m * (1 / (1 - fMRm[2])) - FO2m 
FO2MRm3 = FO2m * (1 / (1 - fMRm[3])) - FO2m
FO2MRm4 = FO2m * (1 / (1 - fMRm[4])) - FO2m
FO2MRm5 = FO2m * (1 / (1 - fMRm[5])) - FO2m

# Total marine biosphere O2 production
FO2m0 = FO2m + FO2MRm0
FO2m1 = FO2m + FO2MRm1
FO2m2 = FO2m + FO2MRm2
FO2m3 = FO2m + FO2MRm3
FO2m4 = FO2m + FO2MRm4
FO2m5 = FO2m + FO2MRm5

#%% Global biosphere

# Total global biosphere O2 production
FO2g0 = FO2t + FO2m0
FO2g1 = FO2t + FO2m1
FO2g2 = FO2t + FO2m2
FO2g3 = FO2t + FO2m3
FO2g4 = FO2t + FO2m4
FO2g5 = FO2t + FO2m5

# Global photorespiration O2 uptake fraction
fPR0 = FO2PRt / FO2g0
fPR1 = FO2PRt / FO2g1
fPR2 = FO2PRt / FO2g2
fPR3 = FO2PRt / FO2g3
fPR4 = FO2PRt / FO2g4
fPR5 = FO2PRt / FO2g5

# Terrestrial Mehler reaction uptake fraction
fMRt0 = FO2MRt / FO2g0
fMRt1 = FO2MRt / FO2g1
fMRt2 = FO2MRt / FO2g2
fMRt3 = FO2MRt / FO2g3
fMRt4 = FO2MRt / FO2g4
fMRt5 = FO2MRt / FO2g5

# Marine Mehler reaction uptake fraction
fMRm0 = FO2MRm0 / FO2g0
fMRm1 = FO2MRm1 / FO2g1
fMRm2 = FO2MRm2 / FO2g2
fMRm3 = FO2MRm3 / FO2g3
fMRm4 = FO2MRm4 / FO2g4
fMRm5 = FO2MRm5 / FO2g5

# Global cytochrome oxidase uptake fraction
fCOX0 = 1 - fPR0 - fMRt0 - fMRm0 
fCOX1 = 1 - fPR1 - fMRt1 - fMRm1 
fCOX2 = 1 - fPR2 - fMRt2 - fMRm2 
fCOX3 = 1 - fPR3 - fMRt3 - fMRm3
fCOX4 = 1 - fPR4 - fMRt4 - fMRm4
fCOX5 = 1 - fPR5 - fMRt5 - fMRm5

# Global alternative oxidase uptake fraction 
fAOX0 = fCOX0 * fAOXCOX
fAOX1 = fCOX1 * fAOXCOX
fAOX2 = fCOX2 * fAOXCOX
fAOX3 = fCOX3 * fAOXCOX
fAOX4 = fCOX4 * fAOXCOX
fAOX5 = fCOX5 * fAOXCOX

# Global cytochrome uptake fraction taking alternative oxidase into account
fCOX0 = fCOX0 - fAOX0
fCOX1 = fCOX1 - fAOX1
fCOX2 = fCOX2 - fAOX2
fCOX3 = fCOX3 - fAOX3
fCOX4 = fCOX4 - fAOX4
fCOX5 = fCOX5 - fAOX5

#%% Calculate tGA

# Relevant thetas

# COX, PR, MR, and AOX Thetas from Hebrew University
tCOXHU = 0.516 # tCOX = 0.516 from Helman2005
tPRHU = 0.512 # tPR = 0.506 from Angert2003 and 0.512 from Helman2005
tMRtHU = 0.526 # tMR = 0.526 for pea from Helman2005
tMRmHU = 0.497 # tMR = 0.497 for Synechocystis from Helman2005
tAOXHU = 0.514 # tAOX = 0.514 from Angert2003

# COX, PR, MR, and AOX Thetas from Rice University
tCOXRU = 0.520 # tCOX = 0.520 from Ash 2019
tPRRU = tPRHU + 0.004 # tPR = 0.516 offset from Helman2005
tMRtRU = tMRtHU + 0.004 # tMR = 0.530 offset from pea from Helman2005
tMRmRU = tMRmHU + 0.004 # tMR = 0.501 offset from Synechocystis from Helman2005
tAOXRU = tAOXHU + 0.004 # tAOX = 0.518 offset from Angert2003

# Function determining tGA
def tGA(fCOX, tCOX, fPR, tPR, fMRt, tMRt, fMRm, tMRm, fAOX, tAOX):
    tGA = fCOX * tCOX + fPR * tPR + fMRt * tMRt + fMRm * tMRm + fAOX * tAOX
    return tGA;

# tGA for Hebrew University thetas

# For 0, 10, 20, 30, 40, and 50% marine Mehler fractions
tGAHU0 = tGA(fCOX0, tCOXHU, fPR0, tPRHU, fMRt0, tMRtHU, fMRm0, tMRmHU,
             fAOX0, tAOXHU)
tGAHU1 = tGA(fCOX1, tCOXHU, fPR1, tPRHU, fMRt1, tMRtHU, fMRm1, tMRmHU,
             fAOX1, tAOXHU)
tGAHU2 = tGA(fCOX2, tCOXHU, fPR2, tPRHU, fMRt2, tMRtHU, fMRm2, tMRmHU,
             fAOX2, tAOXHU)
tGAHU3 = tGA(fCOX3, tCOXHU, fPR3, tPRHU, fMRt3, tMRtHU, fMRm3, tMRmHU,
            fAOX3, tAOXHU)
tGAHU4 = tGA(fCOX4, tCOXHU, fPR4, tPRHU, fMRt4, tMRtHU, fMRm4, tMRmHU,
            fAOX4, tAOXHU)
tGAHU5 = tGA(fCOX5, tCOXHU, fPR5, tPRHU, fMRt5, tMRtHU, fMRm5, tMRmHU,
            fAOX5, tAOXHU)

# TGA for Rice university thetas

# For 0, 10, 20, 30, 40, and 50% marine Mehler fractions
tGARU0 = tGA(fCOX0, tCOXRU, fPR0, tPRRU, fMRt0, tMRtRU, fMRm0, tMRmRU,
             fAOX0, tAOXRU)
tGARU1 = tGA(fCOX1, tCOXRU, fPR1, tPRRU, fMRt1, tMRtRU, fMRm1, tMRmRU,
             fAOX1, tAOXRU)
tGARU2 = tGA(fCOX2, tCOXRU, fPR2, tPRRU, fMRt2, tMRtRU, fMRm2, tMRmRU,
             fAOX2, tAOXRU)
tGARU3 = tGA(fCOX3, tCOXRU, fPR3, tPRRU, fMRt3, tMRtRU, fMRm3, tMRmRU,
             fAOX3, tAOXRU)
tGARU4 = tGA(fCOX4, tCOXRU, fPR4, tPRRU, fMRt4, tMRtRU, fMRm4, tMRmRU,
             fAOX4, tAOXRU)
tGARU5 = tGA(fCOX5, tCOXRU, fPR5, tPRRU, fMRt5, tMRtRU, fMRm5, tMRmRU,
             fAOX5, tAOXRU)

# Append all
tGA = np.append(tGAHU0, tGAHU1)
tGA = np.append(tGA, tGAHU2)
tGA = np.append(tGA, tGAHU3)
tGA = np.append(tGA, tGAHU4)
tGA = np.append(tGA, tGAHU5)
tGA = np.append(tGA, tGARU0)
tGA = np.append(tGA, tGARU1)
tGA = np.append(tGA, tGARU2)
tGA = np.append(tGA, tGARU3)
tGA = np.append(tGA, tGARU4)
tGA = np.append(tGA, tGARU5)

#%% D17O2t as function of fCOX, fPR, and fMR

# Lists for outputs of for loop
sol = []
fracfluxsol = []

for tGAi in tGA:
    
    # Terrestrial fraction of global primary production from Field 1998
    ft = 0.6
        
    # Initial moles of combined isotopologues for relevant species
    Os0 = 1 # O strat – initial moles, any small value works
    O1Ds0 = 1 # O(1D) strat – initial moles, any small value works
    O2s0 = 1 # O2 strat – initial moles, any small value works
    CO2s0 = 4.8e15 # CO2 strat – initial moles, fixes mixing ratio of CO2
    O3s0 = 1 # O3 strat – initial moles, any small value works
    O2t0 = 1 # O2 trop – initial moles, any small value works
    CO2t0 = 4.8e16 # CO2 trop – intial moles (270 ppm), 400ppm = 7.2e16
    O2b0 = 1.83e19 # O2 bio total – from H2O initial moles, not used when H2O is infinite
    
    # Number of oxygens in each species
    isok1 = pd.Series(np.array([1, 1, 2, 2, 3, 2, 2, 2, 2]), 
                      index = np.array(['Os', 'O1Ds', 'O2s', 'CO2s', 'O3s', 'O2t',
                                        'CO2t','O2b', 'O2g']))
    
    # Multiplier for converting atomic ratios to isotopologue ratios for each species
    isok2 = pd.Series(np.array([1, 1, 2, 1, 3, 2, 1, 2, 2]), 
                      index = np.array(['Os', 'O1Ds', 'O2s', 'CO2s', 'O3s', 'O2t',
                                        'CO2t','O2b', 'O2g']))
    
    # SMOW ratios
    rQSMOW = 0.0020052
    rXSMOW = 0.0003799
    
    # Calculate 18O/16O and 17O/16O ratios for each species
    
    # Z is either 18O (Q) or 17O (X)
    def xxZ(isok2, rZ):
        xxZ = (isok2 * rZ) / (isok2 * rQSMOW + isok2 * rXSMOW + 1)
        return xxZ
    
    for species in isok2:
        xxQ = xxZ(isok2, rQSMOW)
        xxX = xxZ(isok2, rXSMOW)
        
    # Calculate initial moles of species with Q and X using isotopologue ratios
    Qs0 = xxQ[0] * Os0 # Q strat - initial moles
    Xs0 = xxX[0] * Os0 # X strat - initial moles
    Q1Ds0 = xxQ[1] * O1Ds0 # Q1D strat - inital moles
    X1Ds0 = xxX[1] * O1Ds0 # X1D strat - intial moles
    OQs0 = xxQ[2] * O2s0 # OQ strat - inital moles
    OXs0 = xxX[2] * O2s0 # OX strat - initial moles
    COQs0 = xxQ[3] * CO2s0 # COQ strat - initial moles
    COXs0 = xxX[3] * CO2s0 # COX strat - initial moles
    OOQs0 = xxQ[4] * O3s0 # OOQ strat - initial moles
    OOXs0 = xxX[4] * O3s0 # OOX strat - initial moles
    OQt0 = xxQ[5] * O2t0 # OQ trop - initial moles
    OXt0 = xxX[5] * O2t0 # OX trop - initial moles
    COQt0 = xxQ[6] * CO2t0 # COQ trop - initial moles
    COXt0 = xxX[6] * CO2t0 # COX trop - initial moles
    OQb0 = xxQ[7] * O2b0 # OQ bio - initial moles
    OXb0 = xxX[7] * O2b0 # OX bio - initial moles
    
    # Initial moles of common isotopologues (i.e. no Q or X)
    Os0 = Os0 - Qs0 - Xs0 # O strat - initial moles
    O1Ds0 = O1Ds0 - Q1Ds0 - X1Ds0 # O1D strat - initial moles
    O2s0 = O2s0 - OQs0 - OXs0 # O2 strat - initial moles
    CO2s0 = CO2s0 - COQs0 - COXs0 # CO2 strat - inital moles
    O3s0 = O3s0 - OOQs0 - OOXs0 # O3 strat - initial moles
    O2t0 = O2t0 - OQt0 - OXt0 # O2 strat - initial moles
    CO2t0 = CO2t0 - COQt0 - COXt0 # CO2 strat - initial moles
    O2b0 = O2b0 - OQb0 - OXb0 # O2 bio - initial moles
    
    #%% Constants relevant to calculations
    
    # Updated fractionation factors
    alphaCO2H2O = 1.041 # 45 oC from Beck et al. 2005
    tequil = 0.528 # nominal TOI equilibration slope
    tGA = tGAi # global average oxygen uptake theta from Young2014
    tevap = 0.5154 # evapotranspiration theta from Landais2006
    tphoto = 0.525 # photosynthetic theta from LuzBarkan2011
    isoevap = 1.005 # isotopic enrichment in water by evapotranspiration from West2008
    isophoto = 1.003389 # net isotopic enrichment of photosynthetic O2 from Luz2014
    alphart = (1 / isoevap) * (1 / 1.01695)  # terrestrial respiration fractionation factor
    alpharm = (1 / isophoto) * (1 / 1.025) # marine respiration fractionation factor
    
    # Moles of air in stratosphere and troposphere
    airs = 1.8e19 # moles of air in the stratosphere
    airt = airs * 10 # moles of air in the troposphere is 10x that of the stratosphere
    
    # Volume of stratosphere and troposphere in cm^3
    vs = 2.8e25 # volume of stratosphere
    vt = vs * 10 # volume of troposphere is 10x that of the stratosphere
    
    # Constants for calculating reaction rates
    secyear = 31536000 # number of seconds in a year
    avo = 6.0221409e23 # Avogadro's number
    
    #%% Calculate reaction rates
    
    # Equations used to calculate reaction rates
    
    # Convert reaction rate from units of cm^3/s to 1 / (mol yr)
    def tomol(cm3s):
        tomol = cm3s * secyear * (avo / vs)
        return tomol
    
    # Low pressure limit rate constant calculation from JPL 19-5
    def lpk(ko298, T, n):
        k = ko298 * (T / 298) ** (-1 * n)
        return k;
    
    # Arrhenius rate constant expression from JPL 19-5
    def Arr(A, T, ER):
        k = A * np.exp(-1 * ER * (1 / T))
        return k;
    
    # O2 + O isotope exchange expression from Fleurat-Lessard2003
    def O2O(T):
        O2O = 2.7e-12 * (300 / T) ** 0.9
        return O2O;
    
    # Rate constants for transport between boxes
    # 1 - troposphere, 2 - biosphere/hydrosphere, 3 - geosphere, 4 - stratosphere
    # k12 represents 1 (trop) -> 2 (bio/hydro). similar notation for others
    GPPx = 1
    k12 = 0.0008 * GPPx # respiration rate constant yr^-1
    k12t = ft * k12 # respiration rate constant from bio terrestrial yr^-1
    k12m = (1 - ft) * k12 # respiration rate constant from bio marine yr^-1
    k21 = 0.001653 * GPPx # photosynthesis rate constant yr^-1
    k21t = ft * k21 # photosynthesis rate constant from bio terrestrial yr^-1
    k21m = (1 - ft) * k21 # photosynthesis rate constant from bio marine yr^-1
    k41 = 1 # stat-trop mixing rate constant yr^-1
    k14 = k41 / 10 # trop-strat mixing rate constant yr^-1
    
    #%% Calculate reaction rates using JPL 19-5 data without reduced mass scaling
    KMIF = 1.065 # MIF for O3 formation
    K1 = 1.109e-12 * secyear * 1 # O2 + PHO -> O + O 1/(yr mol)
    K2 = 1.109e-12 * secyear * 1 # OQ + PHO -> Q + O 1/(yr mol)
    K3 = 1.109e-12 * secyear * 1 # OX + PHO -> X + O 1/(yr mol)
    k4i = lpk(6.1e-34, 220, 2.4) * 8.3e17
    K4 = tomol(k4i) * 1 # O2 + O -> O3 1/(yr mol)
    K5 = KMIF * tomol(k4i) * 1 # O2 + Q -> OOQ 1/(yr mol)
    K6 = KMIF * tomol(k4i) * 1 # O2 + X -> OOX 1/(yr mol)
    K7 = KMIF * tomol(k4i) * 1 # OQ + O -> OOQ 1/(yr mol)
    K8 = KMIF * tomol(k4i) * 1 # OX + O -> OOX 1/(yr mol)
    K9 = 2.96e-4 * secyear * 1 # O3 + PHO -> O2 + O 1/(yr mol)
    K10 = 2.96e-4 * secyear * (1/3) * 1 # OOQ + PHO -> O2 + Q 1/(yr mol)
    K11 = 2.96e-4 * secyear * (2/3) * 1 # OOQ + PHO -> OQ + O 1/(yr mol)
    K12 = 2.96e-4 * secyear * (1/3) * 1 # OOX + PHO -> O2 + X 1/(yr mol)
    K13 = 2.96e-4 * secyear * (2/3) * 1 # OOX + PHO -> OX + O 1/(yr mol)
    K14 = 5.01e-4 * secyear * 1 # O3 + PHO -> O2 + O1D 1/(yr mol)
    K15 = 5.01e-4 * secyear * (1/3) * 1 # OOQ + PHO -> O2 + Q1D 1/(yr mol)
    K16 = 5.01e-4 * secyear * (2/3) * 1 # OOQ + PHO -> OQ + O1D 1/(yr mol)
    K17 = 5.01e-4 * secyear * (1/3) * 1 # OOX + PHO -> O2 + X1D 1/(yr mol)
    K18 = 5.01e-4 * secyear * (2/3) * 1 # OOX + PHO -> OX + O1D 1/(yr mol)
    k19i = Arr(8e-12, 220, 2060)
    K19 = tomol(k19i) * 1 # O3 + O -> O2 + O2 mol/yr
    K20 = tomol(k19i) * 1 # OOQ + O -> O2 + OQ 1/(yr mol)
    K21 = tomol(k19i) * 1  # OOX + O -> O2 + OX 1/(yr mol)
    K22 = tomol(k19i) * 1 # O3 + Q -> O2 + OQ 1/(yr mol)
    K23 = tomol(k19i) * 1  # O3 + X -> O2 + OX 1/(yr mol)
    k24i = Arr(2.4e-10, 220, 0)
    K24 = tomol(k24i) * 1  # O3 + O1D -> O2 + O2 1/(yr mol)
    K25 = tomol(k24i) * (1/2) * 1 # OOQ + O1D -> O2 + OQ 1/(yr mol)
    K26 = tomol(k24i) * (1/2) * 1 # OOX + O1D -> O2 + OX 1/(yr mol)
    K27 = tomol(k24i) * (1/2) * 1 # O3 + Q1D -> O2 + OQ 1/(yr mol)
    K28 = tomol(k24i) * (1/2) * 1 # O3 + X1D -> O2 + OX 1/(yr mol)
    K29 = tomol(k24i) * (1/2) * 1 # O3 + O1D -> O2 + O + O 1/(yr mol)
    K30 = tomol(k24i) * (1/2) * (1/2) * 1 # OOQ + O1D -> O2 + O + Q 1/(yr mol)
    K31 = tomol(k24i) * (1/2) * (1/2) * 1 # OOQ + O1D -> OQ + O + O 1/(yr mol)
    K32 = tomol(k24i) * (1/2) * (1/2) * 1 # OOX + O1D -> O2 + O + X 1/(yr mol)
    K33 = tomol(k24i) * (1/2) * (1/2) * 1 # OOX + O1D -> OX + O + O 1/(yr mol)
    K34 = tomol(k24i) * (1/2) * (1/2) * 1 # O3 + Q1D -> O2 + O + Q 1/(yr mol)
    K35 = tomol(k24i) * (1/2) * (1/2) * 1 # O3 + Q1D -> OQ + O + O 1/(yr mol)
    K36 = tomol(k24i) * (1/2) * (1/2) * 1 # O3 + X1D -> O2 + O + X 1/(yr mol)
    K37 = tomol(k24i) * (1/2) * (1/2) * 1 # O3 + X1D -> OX + O + O 1/(yr mol)
    k38i = Arr(3.3e-11, 220, -55)
    K38 = tomol(k38i) * ((0.78 + 0.21) / 0.21) * 1 # O2 + O1D -> O2 + O mol/yr
    K39 = tomol(k38i) * ((0.78 + 0.21) / 0.21) * 1 # O2 + Q1D -> O2 + Q 1/(yr mol)
    K40 = 0 # O2 + Q1D -> OQ + O 1/(yr mol)
    K41 = tomol(k38i) * ((0.78 + 0.21) / 0.21) * 1 # O2 + X1D -> O2 + X 1/(yr mol)
    K42 = 0 # O2 + X1D -> OX + O 1/(yr mol)
    k43i = O2O(220)
    K43 = tomol(k43i) * 1 # O2 + Q -> OQ + O 1/(yr mol)
    K44 = tomol(k43i) * (1/2) * 1 # OQ + O -> O2 + Q 1/(yr mol)
    K45 = tomol(k43i) * 1 # O2 + X -> OX + O 1/(yr mol)
    K46 = tomol(k43i) * (1/2) * 1 # OX + O -> O2 + X 1/(yr mol)
    k47i = Arr(7.5e-11, 220, -115)
    TNF = 0.027
    K47 = tomol(k47i) * TNF # CO2 + O1D -> CO2 + O 1/(yr mol)
    K48 = tomol(k47i) * (1/3) * TNF # CO2 + Q1D -> CO2 + Q 1/(yr mol)
    K49 = tomol(k47i) * (2/3) * TNF # CO2 + Q1D -> COQ + O 1/(yr mol)
    K50 = tomol(k47i) * (1/3) * TNF # COQ + O1D -> CO2 + Q 1/(yr mol)
    K51 = tomol(k47i) * (2/3) * TNF # COQ + O1D -> COQ + O 1/(yr mol)
    K52 = tomol(k47i) * (1/3) * TNF # CO2 + X1D -> CO2 + X 1/(yr mol)
    K53 = tomol(k47i) * (2/3) * TNF # CO2 + X1D -> COX + O 1/(yr mol)
    K54 = tomol(k47i) * (1/3) * TNF # COX + O1D -> CO2 + X 1/(yr mol)
    K55 = tomol(k47i) * (2/3) * TNF # COX + O1D -> COX + O 1/(yr mol)
    
    # Calculated hydrosphere rate constants
    
    # Turnover times for terrestrial and marine biospheres
    ttT = 1 # years
    ttM = 1 # Years
    kr10t = 1 / ttT * 1 # COQ + H2O -> CO2 + H2Q 1/(yr mol) terrestrial
    kr10m = 1 / ttM * 1 # COQ + H2O -> CO2 + H2Q 1/(yr mol) marine
    kr9t = kr10t * alphaCO2H2O * isoevap * rQSMOW * 1 # CO2 + H2Q -> COQ + H2O 1/(yr mol) terrestrial
    kr9m = kr10m * alphaCO2H2O * rQSMOW * 1 # CO2 + H2Q -> COQ + H2O 1/(yr mol) marine
    kr12t = kr10t * 1 # COX + H2O -> CO2 + H2X 1/(yr mol) terrestrial 
    kr12m = kr10m * 1 # COX + H2O -> CO2 + H2X 1/(yr mol)
    kr11t = kr12t * alphaCO2H2O ** tequil * isoevap ** tevap * rXSMOW * 1 # CO2 + H2X -> COX + H2O 1/(yr mol) terrestrial
    kr11m = kr12m * alphaCO2H2O ** tequil * rXSMOW * 1 # COX + H2O -> COX + H2X 1/(yr mol) marine
    
    #%% Solve the system of differential equations
    
    # Function containing all the ODEs
    def f(y, t):
        Osi = y[0]
        Xsi = y[1]
        O1Dsi = y[2]
        Qsi = y[3] 
        X1Dsi = y[4]
        Q1Dsi = y[5]
        O2si = y[6]
        OXsi = y[7]
        OQsi = y[8]
        CO2si = y[9]
        COXsi = y[10]
        COQsi = y[11]
        O3si = y[12]
        OOXsi = y[13]
        OOQsi = y[14]
        O2ti = y[15]
        OQti = y[16]
        OXti = y[17]
        CO2ti = y[18]
        COQti = y[19]
        COXti = y[20]
        O2bi = y[21]
        OQbi = y[22]
        OXbi = y[23]
       
        # k(species)I represents input flux, k(species)O represents output flux
        
        # Stratosphere ODEs
        
        # O Stratosphere
        kOsI = (2 * K1 * O2si + K2 * OQsi + K3 * OXsi + K9 * O3si + K11 * OOQsi +
                K13 * OOXsi + 2 * K29 * O3si * O1Dsi + K30 * OOQsi * O1Dsi +
                2 * K31 * OOQsi * O1Dsi + K32 * OOXsi * O1Dsi +
                2 * K33 * OOXsi * O1Dsi + K34 * O3si * Q1Dsi +
                2 * K35 * O3si * Q1Dsi + K36 * O3si * X1Dsi +
                2 * K37 * O3si * X1Dsi + K38 * O2si * O1Dsi +
                K40 * O2si * Q1Dsi + K42 * O2si * X1Dsi + K43 * O2si * Qsi +
                K45 * O2si * Xsi + K47 * CO2si * O1Dsi + K49 * CO2si * Q1Dsi +
                K51 * COQsi * O1Dsi + K53 * CO2si * X1Dsi + K55 * COXsi * O1Dsi)
        kOsO = (K4 * O2si + K7 * OQsi + K8 * OXsi + K19 * O3si + K20 * OOQsi +
                K21 * OOXsi + K44 * OQsi + K46 * OXsi)
        dOs = kOsI - Osi * kOsO
        #dOs = 0
        
        # X (17O) Stratosphere 
        kXsI = (K3 * OXsi + K12 * OOXsi + K32 * OOXsi * O1Dsi +
                K36 * O3si * X1Dsi + K41 * O2si * X1Dsi + K46 * OXsi * Osi +
                K52 * CO2si * X1Dsi + K54 * COXsi * O1Dsi)
        kXsO = (K6 * O2si + K23 * O3si + K45 * O2si)
        dXs = kXsI - Xsi * kXsO
        #dXs = 0
        
        # Q (18O) Stratosphere
        kQsI = (K2 * OQsi + K10 * OOQsi + K30 * OOQsi * O1Dsi +
                K34 *  O3si * Q1Dsi + K39 * O2si * Q1Dsi + K44 * OQsi * Osi +
                K48 * CO2si * Q1Dsi + K50 * COQsi * O1Dsi)
        kQsO = (K5 * O2si + K22 * O3si + K43 * O2si)
        dQs = kQsI - Qsi * kQsO
        #dQs = 0
        
        # O(1D) Stratosphere
        kO1DsI = (K14 * O3si + K16 * OOQsi + K18 * OOXsi)
        kO1DsO = (K24 * O3si + K25 * OOQsi + K26 * OOXsi + K29 * O3si +
                  K30 * OOQsi + K31 * OOQsi + K32 * OOXsi + K33 * OOXsi +
                  K38 * O2si + K47 * CO2si + K50 * COQsi + K51 * COQsi +
                  K54 * COXsi + K55 * COXsi)
        dO1Ds = kO1DsI - O1Dsi * kO1DsO 
        #dO1Ds = 0
        
        # X(1D) (17O(1D)) Stratosphere
        kX1DsI = (K17 * OOXsi)
        kX1DsO = (K28 * O3si + K36 * O3si + K37 * O3si + K41 * O2si +
                  K42 * O2si + K52 * CO2si + K53 * CO2si)
        dX1Ds = kX1DsI - X1Dsi * kX1DsO
        #dX1Ds = 0
        
        # Q(1D) (18O(1D)) Stratosphere
        kQ1DsI = (K15 * OOQsi)
        kQ1DsO = (K27 * O3si + K34 * O3si + K35 * O3si + K39 * O2si +
                  K40 * O2si + K48 * CO2si + K49 * CO2si)
        dQ1Ds = kQ1DsI - Q1Dsi * kQ1DsO
        #dQ1Ds = 0
        
        # O2 Stratosphere
        kO2sI = (K9 * O3si + K10 * OOQsi + K12 * OOXsi + K14 * O3si +
                 K15 * OOQsi + K17 * OOXsi + 2 * K19 * O3si * Osi +
                 K20 * OOQsi * Osi + K21 * OOXsi * Osi + K22 * O3si * Qsi +
                 K23 * O3si * Xsi + 2 * K24 * O3si * O1Dsi +
                 K25 * OOQsi * O1Dsi + K26 * OOXsi * O1Dsi +
                 K27 * O3si * Q1Dsi + K28 * O3si * X1Dsi + K29 * O3si * O1Dsi +
                 K30 * OOQsi * O1Dsi + K32 * OOXsi * O1Dsi + K34 * O3si * Q1Dsi +
                 K36 * O3si * X1Dsi + K38 * O2si * O1Dsi + K39 * O2si * Q1Dsi +
                 K41 * O2si * X1Dsi + K44 * OQsi * Osi + K46 * OXsi * Osi +
                 k14 * O2ti)
        kO2sO = (K1 + K4 * Osi + K5 * Qsi + K6 * Xsi + K38 * O1Dsi +
                 K39 * Q1Dsi + K40 * Q1Dsi + K41 * X1Dsi + K42 * X1Dsi +
                 K43 * Qsi + K45 * Xsi + k41)
        dO2s = kO2sI - O2si * kO2sO
        #dO2s = 0
        
        # OX (O17O) Stratosphere
        kOXsI = (K13 * OOXsi + K18 * OOXsi + K21 * OOXsi * Osi +
                 K23 * O3si * Xsi + K26 * OOXsi * O1Dsi + K28 * O3si * X1Dsi +
                 K33 * OOXsi * O1Dsi + K37 * O3si * X1Dsi + K42 * O2si * X1Dsi +
                 K45 * O2si * Xsi + k14 * OXti)
        kOXsO = (K3 + K8 * Osi + K46 * Osi + k41)
        dOXs = kOXsI - OXsi * kOXsO
        #dOXs = 0
        
        # OQ (O18O) Stratosphere
        kOQsI = (K11 * OOQsi + K16 * OOQsi + K20 * OOQsi * Osi +
                 K22 * O3si * Qsi + K25 * OOQsi * O1Dsi + K27 * O3si * Q1Dsi +
                 K31 * OOQsi * O1Dsi + K35 *  O3si * Q1Dsi + K40 * O2si * Q1Dsi +
                 K43 * O2si * Qsi + k14 * OQti)
        kOQsO = (K2 + K7 * Osi + K44 * Osi + k41)
        dOQs = kOQsI - OQsi * kOQsO
        #dOQs = 0
        
        # CO2 Stratosphere
        kCO2sI = (K47 * CO2si * O1Dsi + K48 * CO2si * Q1Dsi + K50 * COQsi * O1Dsi +
                  K52 * CO2si * X1Dsi + K54 * COXsi * O1Dsi + k14 * CO2ti)
        kCO2sO = (K47 * O1Dsi + K48 * Q1Dsi + K49 * Q1Dsi + K52 * X1Dsi +
                  K53 * X1Dsi + k41)
        dCO2s = kCO2sI - CO2si * kCO2sO
        #dCO2s = 0
        
        # COX (CO17O) Stratosphere
        kCOXsI = (K53 * CO2si * X1Dsi + K55 * COXsi * O1Dsi + k14 * COXti)
        kCOXsO = (K54 * O1Dsi + K55 * O1Dsi + k41)
        dCOXs = kCOXsI - COXsi * kCOXsO
        #dCOXs = 0
        
        # COQ (CO18O) Stratosphere
        kCOQsI = (K49 * CO2si * Q1Dsi + K51 * COQsi * O1Dsi + k14 * COQti)
        kCOQsO = (K50 * O1Dsi + K51 * O1Dsi + k41)
        dCOQs = kCOQsI - COQsi * kCOQsO
        #dCOQs = 0
        
        # O3 Stratosphere
        kO3sI = (K4 * O2si * Osi)
        kO3sO = (K9 + K14 + K19 * Osi + K22 * Qsi + K23 * Xsi +
                 K24 * O1Dsi + K27 * Q1Dsi + K28 * X1Dsi + K29 * O1Dsi +
                 K34 * Q1Dsi + K35 * Q1Dsi + K36 * X1Dsi + K37 * X1Dsi)
        dO3s = kO3sI - O3si * kO3sO
        #dO3s = 0
        
        # OOX (OO17O) Stratosphere
        kOOXsI = (K6 * O2si * Xsi + K8 * OXsi * Osi)
        kOOXsO = (K12 + K13 + K17 + K18 + K21 * Osi + K26 * O1Dsi +
                  K32 * O1Dsi + K33 * O1Dsi)
        dOOXs = kOOXsI - OOXsi * kOOXsO
        #dOOXs = 0
        
        # OOQ (OO18O) Stratosphere
        kOOQsI = (K5 * O2si * Qsi + K7 * OQsi * Osi)
        kOOQsO = (K10 + K11 + K15 + K16 + K20 * Osi + K25 * O1Dsi +
                  K30 * O1Dsi + K31 * O1Dsi)
        dOOQs = kOOQsI - OOQsi * kOOQsO
        #dOOQs = 0
        
        # Troposphere ODEs
        
        # O2 Troposphere
        kO2tI = (k41 * O2si + k21t * O2bi + k21m * O2bi)
        kO2tO = (k12t + k12m + k14)
        dO2t = kO2tI - O2ti * kO2tO
        #dO2t = 0
        
        # OX (O17O) Troposphere
        kOXtI = (k41 * OXsi + k21t * OXbi + k21m * OXbi)
        kOXtO = (k12t * alphart ** tGA + k12m * alpharm ** tGA + k14)
        dOXt = kOXtI - OXti * kOXtO
        #dOXt = 0
        
        # OQ (O18O) Troposphere
        kOQtI = (k41 * OQsi + k21t * OQbi + k21m * OQbi)
        kOQtO = (k12t * alphart + k12m * alpharm + k14)
        dOQt = kOQtI - OQti * kOQtO
        #dOQt = 0
        
        # CO2 Troposphere
        kCO2tI = (kr10t * COQti + kr10m * COQti + kr12t * COXti + kr12m * COXti + k41 * CO2si)
        kCO2tO = (kr9t + kr9m + kr11t + kr11m + k14)
        dCO2t = kCO2tI - CO2ti * kCO2tO
        #dCO2t = 0
        
        # COX (CO17O) Troposphere
        kCOXtI = (kr11t * CO2ti + kr11m * CO2ti + k41 * COXsi) 
        kCOXtO = (kr12t + kr12m + k14)
        dCOXt = kCOXtI - COXti * kCOXtO
        #dCOXt = 0
        
        # COQ (CO18O) Troposphere
        kCOQtI = (kr9t * CO2ti + kr9m * CO2ti  + k41 * COQsi)
        kCOQtO = (kr10t + kr10m + k14) 
        dCOQt = kCOQtI - COQti * kCOQtO
        #dCOQt = 0
        
        # Biosphere ODEs (equations not used because we assume an infinite reservoir)
        
        # O2 Biosphere
        kO2bI = (k12t * O2ti + k12m * O2ti)
        kO2bO = (k21t + k21m)
        dO2b = kO2bI - O2bi *  kO2bO
        dO2b = 0
        
        # OX (17OO) Biosphere
        kOXbI = (k12t * alphart ** tGA * OXti + k12m * alpharm **tGA * OXti)
        kOXbO = (k21t + k21m)
        dOXb = kOXbI - OXbi * kOXbO
        dOXb = 0
        
        # OQ (18OO) Biosphere
        kOQbI = (k12t * alphart * OQti + k12m * alpharm * OQti)
        kOQbO = (k21t + k21m)
        dOQb = kOQbI - OQbi * kOQbO
        dOQb = 0
        
        return np.array([dOs, dXs, dO1Ds, dQs, dX1Ds, dQ1Ds, dO2s, dOXs, dOQs,
                         dCO2s, dCOXs, dCOQs, dO3s, dOOXs, dOOQs, dO2t, dOQt,
                         dOXt, dCO2t, dCOQt, dCOXt,dO2b, dOQb, dOXb])
    
    # Initial conditions for solver
    y0 = np.array([Os0, Xs0, O1Ds0, Qs0, X1Ds0, Q1Ds0, O2s0, OXs0, OQs0, CO2s0,
                   COXs0, COQs0, O3s0, OOXs0, OOQs0, O2t0, OQt0, OXt0, CO2t0,
                   COQt0, COXt0, O2b0, OQb0, OXb0])
    
    
    # Time grid
    t = np.arange(0, 10e5, 0.1)
    
    # Order of species in molar output
    moleso = np.array(['Os', 'Xs', 'O1Ds', 'Qs', 'X1Ds', 'Q1Ds', 'O2s', 'OXs',
                       'OQs', 'CO2s', 'COXs', 'COQs', 'O3s', 'OOXs', 'OOQs',
                       'O2t', 'OQt', 'OXt', 'CO2t', 'COQt', 'COXt', 'O2b',
                       'OQb', 'OXb'])
    
    # Solve the DEs
    moles = odeint(f, y0, t, mxstep = 1000000)
        
    # Set up output dataframe
    moles = pd.DataFrame(moles, columns = moleso)
    moles = pd.concat([moles.tail(1)]).reset_index(drop=True)
    moles = moles.transpose().reset_index(drop=True)
    moles['Index'] = moleso
    moles = moles.set_index('Index')
    
    #%% Isotopes output
    
    # Functions that calculate atomic isotopologue ratios (18O/16O and 17O/16O)
    # Z is either 18O (Q) or 17O (X) and O refers to non rare isotopologue
    def atomZ(isok2, Z, O):
        atomZ = (1 / isok2) * (Z / O)
        return atomZ
    
    # Use atomic isotopologue ratios to calculate delta values relative to SMOW
    def deltaZ(atomZ, rZ):
        deltaZ = 1000 * np.log(atomZ / rZ)
        return deltaZ
    
    # Use delta values to calculate cap17
    def capD(deltaX, deltaQ):
        capD = deltaX - tequil * deltaQ
        return capD
    
    # Calculate the 18O/16O and 17O/16O ratios of molecules at end of model run
    R18_Os = atomZ(isok2[0], moles.loc['Qs'], moles.loc['Os'])
    R17_Os = atomZ(isok2[0], moles.loc['Xs'], moles.loc['Os'])
    R18_O1Ds = atomZ(isok2[1], moles.loc['Q1Ds'], moles.loc['O1Ds'])
    R17_O1Ds = atomZ(isok2[1], moles.loc['X1Ds'], moles.loc['O1Ds'])
    R18_O2s = atomZ(isok2[2], moles.loc['OQs'], moles.loc['O2s'])
    R17_O2s = atomZ(isok2[2], moles.loc['OXs'], moles.loc['O2s'])
    R18_CO2s = atomZ(isok2[3], moles.loc['COQs'], moles.loc['CO2s'])
    R17_CO2s = atomZ(isok2[3], moles.loc['COXs'], moles.loc['CO2s'])
    R18_O3s = atomZ(isok2[4], moles.loc['OOQs'], moles.loc['O3s'])
    R17_O3s = atomZ(isok2[4], moles.loc['OOXs'], moles.loc['O3s'])
    R18_O2t = atomZ(isok2[5], moles.loc['OQt'], moles.loc['O2t'])
    R17_O2t = atomZ(isok2[5], moles.loc['OXt'], moles.loc['O2t'])
    R18_CO2t = atomZ(isok2[6], moles.loc['COQt'], moles.loc['CO2t'])
    R17_CO2t = atomZ(isok2[6], moles.loc['COXt'], moles.loc['CO2t'])
    
    # Calculate the d17O and d18O of species.
    d18_Os = deltaZ(R18_Os, rQSMOW)
    d17_Os = deltaZ(R17_Os, rXSMOW)
    d18_O1Ds = deltaZ(R18_O1Ds, rQSMOW)
    d17_O1Ds = deltaZ(R17_O1Ds, rXSMOW)
    d18_O2s = deltaZ(R18_O2s, rQSMOW)
    d17_O2s = deltaZ(R17_O2s, rXSMOW)
    d18_CO2s = deltaZ(R18_CO2s, rQSMOW)
    d17_CO2s = deltaZ(R17_CO2s, rXSMOW)
    d18_O3s = deltaZ(R18_O3s, rQSMOW)
    d17_O3s = deltaZ(R17_O3s, rXSMOW)
    d18_O2t = deltaZ(R18_O2t, rQSMOW)
    d17_O2t = deltaZ(R17_O2t, rXSMOW)
    d18_CO2t = deltaZ(R18_CO2t, rQSMOW)
    d17_CO2t = deltaZ(R17_CO2t, rXSMOW)
    
    # Calculate D17 at end of model run
    D17_Os = capD(d17_Os, d18_Os)
    D17_O1Ds = capD(d17_O1Ds, d18_O1Ds)
    D17_O2s = capD(d17_O2s, d18_O2s)
    D17_CO2s = capD(d17_CO2s, d18_CO2s)
    D17_O3s = capD(d17_O3s, d18_O3s)
    D17_O2t = capD(d17_O2t, d18_O2t)
    D17_CO2t = capD(d17_CO2t, d18_CO2t)
    
    # Order of species in ratio output dataframe
    ratioso = np.array(['R18_Os', 'R18_O1Ds', 'R18_O2s', 'R18_CO2s', 'R18_O3s',
                      'R18_O2t', 'R18_CO2t', 'R17_Os', 'R17_O1Ds', 'R17_O2s',
                      'R17_CO2s', 'R17_O3s', 'R17_O2t', 'R17_CO2t'])
    
    # Order of species in isotope output dataframe
    isotopeso = np.array(['d18_Os', 'd18_O1Ds', 'd18_O2s', 'd18_CO2s', 'd18_O3s',
                          'd18_O2t', 'd18_CO2t', 'd17_Os', 'd17_O1Ds', 'd17_O2s',
                          'd17_CO2s', 'd17_O3s', 'd17_O2t', 'd17_CO2t', 'D17_Os',
                          'D17_O1Ds', 'D17_O2s', 'D17_CO2s', 'D17_O3s', 'D17_O2t',
                          'D17_CO2t'])
    
    #%% Mole fraction and flux outputs
    
    # Mole fraction of O2 in troposphere
    def xO2(O2, OX, OQ, air):
        xO2 = (O2 + OX + OQ) / air
        return xO2;
    
    # Mole fraction of CO2 in troposphere and stratosphere
    def xCO2(CO2, COX, COQ, air):
        xCO2 = (CO2 + COX + COQ) / air
        return xCO2;
    
    # Mole fraction of O3 in stratosphere
    def xO3(O3, OOX, OOQ, air):
        xO3 = (O3 + OOX + OOQ) / air
        return xO3;
    
    # Mole fraction of O(1D) in stratosphere
    def xO1D(O1D, X1D, Q1D, air):
        xO1D = (O1D + X1D + Q1D) / air
        return xO1D;
    
    # Mole fraction of O in stratosphere
    def xO(O, X, Q, air):
        xO = (O + X + Q) / air
        return xO;
    
    # D17O XO3 molar flux in per mil moles CO2/yr
    def xJ(xCO2s, capDCO2s, airs, k41):
        xJ = (xCO2s * capDCO2s * airs) / (1 / k41)
        return xJ;
    
    # Calculate mole fraction of O2 in the troposphere
    xO2 = xO2(moles.loc['O2t'], moles.loc['OXt'], moles.loc['OQt'], airt)
    
    # Calculate mole fraction of CO2 in the troposphere
    xCO2t = xCO2(moles.loc['CO2t'], moles.loc['COXt'], moles.loc['COQt'], airt)
    
    # Calculate mole fraction of CO2 in the stratosphere
    xCO2s = xCO2(moles.loc['CO2s'], moles.loc['COXs'], moles.loc['COQs'], airs)
    
    # Calculate mole fraction of O3 in the stratosphere
    xO3s = xO3(moles.loc['O3s'], moles.loc['OOXs'], moles.loc['OOQs'], airs)
    
    # Calculate mole fraction of O(1D) in the stratosphere
    xO1Ds = xO1D(moles.loc['O1Ds'], moles.loc['X1Ds'], moles.loc['Q1Ds'], airs)
    
    # Calculate mole fraction of O in the stratosphere
    xOs = xO(moles.loc['Os'], moles.loc['Xs'], moles.loc['Qs'], airs)
    
    # Calculate D17O flux
    xJ = xJ(xCO2s, D17_CO2s, airs, k41)
    
    # Order of species in mole fraction and isotope flux output
    fracfluxo = np.array(['xO2', 'xCO2t', 'xCO2s', 'xO3s', 'xJ'])
    
    #%% Output dataframes to spreadsheet
    
    # Create empty dataframes for isotope and fraction/flux outputs
    isotopes = pd.DataFrame(np.full((1, 21), 0, dtype=float), columns = isotopeso)
    fracflux = pd.DataFrame(np.full((1, 5), 0, dtype=float), columns = fracfluxo)
    
    # Assign calculated isotope values to dataframe
    isotopes = pd.DataFrame(np.array([d18_Os, d18_O1Ds, d18_O2s, d18_CO2s,
                                      d18_O3s, d18_O2t, d18_CO2t, d17_Os,
                                      d17_O1Ds, d17_O2s, d17_CO2s, d17_O3s,
                                      d17_O2t, d17_CO2t, D17_Os, D17_O1Ds,
                                      D17_O2s, D17_CO2s, D17_O3s, D17_O2t,
                                      D17_CO2t]),
                            index = isotopeso)
    
    # Assign calculated ratio values to dataframe
    ratios = pd.DataFrame(np.array([R18_Os, R18_O1Ds, R18_O2s, R18_CO2s, R18_O3s,
                                    R18_O2t, R18_CO2t, R17_Os, R17_O1Ds, R17_O2s,
                                    R17_CO2s, R17_O3s, R17_O2t, R17_CO2t]),
                          index = ratioso)
    
    # Assign calculated mole fraction and flux values to dataframe
    fracflux = pd.DataFrame(np.array([xO2, xCO2t, xCO2s, xO3s, xJ]),
                            index =fracfluxo)
    
    # Add isotopes and fracflux outputs
    sol.append(isotopes)
    fracfluxsol.append(fracflux)

#%% Plots

## Plotting D17O2t for Hebrew University and Rice University

# Setting up figure parameters
fig1= plt.figure(figsize = (5, 5))
fig1 = fig1.add_subplot(1, 1, 1)
fig1.set(xlim = (0.509, 0.521), ylim = (-0.57, -0.29))
fig1.set_facecolor('#F7F7F7')
#fig1.grid()

# D17O2t for 30, 40, and 50% marine Mehler
D17_tGAHU = []
for i in sol[0:6]:
    D17 = i.loc['D17_O2t'].values
    D17_tGAHU.append(D17)
D17_tGAHU = np.hstack(D17_tGAHU)

D17_tGARU = []
for i in sol[6:12]:
    D17 = i.loc['D17_O2t'].values
    D17_tGARU.append(D17)
D17_tGARU = np.hstack(D17_tGARU)

# d17O and d18O of atmospheric O2 from BarkanLuz2011
d17_O2tLB = -11.883
d18_O2tLB = -23.324

# R17 and R18 of LuzBarkan2011 data
R17LB = rXSMOW / ((d17_O2tLB / 1e3) + 1)
R18LB = rQSMOW / ((d18_O2tLB / 1e3) + 1)

# d'17O and d'18O of LuzBarkan2011 data
d17_O2tLB = deltaZ(R17LB, rXSMOW)
d18_O2tLB = deltaZ(R18LB, rQSMOW)

# D'17O of LuzBarkan2011 data
D17_O2tLB = capD(d17_O2tLB, d18_O2tLB)

# Plotting D17 of Wostbrock 2020
xWB = np.linspace(0.509, 0.521, 11)
yWB = np.array([-0.441] * 11)
fig1.plot(xWB, yWB, color='#77aadd', linestyle='solid', zorder=3.5)
fig1.text(0.5095, -0.441 + 0.007, 'Wostbrock2020', color = '#77aadd', fontsize=8)

# Plotting D17 of LuzBarkan2011
xLB = np.linspace(0.509, 0.521, 11)
yLB = np.array([-0.507] * 11)
fig1.plot(xLB, yLB, color='#ee8866', linestyle='solid', zorder=3.5)
fig1.text(0.51785, -0.507 + 0.007, 'LuzBarkan2011', color = '#ee8866', fontsize=8)

# Plotting Rice University D17O2t as a function of tGA
fig1.scatter(tGARU0, D17_tGARU[0], color='#98CAE1', marker='^',
             zorder=4.5) # 0 % MRm 
fig1.text(tGARU0 - 0.0011, D17_tGARU[0] + 0.01, '0%', fontsize=8,
          color='#98CAE1', zorder=4.5)
fig1.scatter(tGARU1, D17_tGARU[1], color='#84B1D3', marker='^',
             zorder=4.5) # 10 % MRm 
fig1.text(tGARU1 - 0.0011, D17_tGARU[1] + 0.01, '10%', fontsize=8,
          color='#84B1D3', zorder=4.5)
fig1.scatter(tGARU2, D17_tGARU[2], color='#7197C5', marker='^',
             zorder=4.5) # 20 % MRm 
fig1.text(tGARU2 - 0.0011, D17_tGARU[2] + 0.01, '20%', fontsize=8,
          color='#7197C5', zorder=4.5)
fig1.scatter(tGARU3, D17_tGARU[3], color='#5D7EB6', marker='^',
             zorder=4.5) # 30 % MRm 
fig1.text(tGARU3 - 0.0011, D17_tGARU[3] + 0.01, '30%', fontsize=8,
          color='#5D7EB6', zorder=4.5)
fig1.scatter(tGARU4, D17_tGARU[4], color='#4A64A8', marker='^',
             zorder=4.5) # 40 % MRm
fig1.text(tGARU4 - 0.0011, D17_tGARU[4] + 0.01, '40%', fontsize=8,
          color='#4A64A8', zorder=4.5)
fig1.scatter(tGARU5, D17_tGARU[5], color='#364B9A', marker='^',
             zorder=4.5, label='Rice University') # 50 % MRm
fig1.text(tGARU5 - 0.0011, D17_tGARU[5] + 0.01, '50%', fontsize=8,
          color='#364B9A', zorder=4.5)

# Plotting Hebrew University D17O2t as a function of tGA
fig1.scatter(tGAHU0, D17_tGAHU[0], color='#FDB366', marker='o',
             zorder=4.5) # 0 % MRm
fig1.text(tGAHU0 + 0.0004, D17_tGAHU[0] - 0.015, '0%', fontsize=8,
          color='#FDB366', zorder=4.5)
fig1.scatter(tGAHU1, D17_tGAHU[1], color='#EB8F59', marker='o',
             zorder=4.5) # 10 % MRm
fig1.text(tGAHU1 + 0.0004, D17_tGAHU[1] - 0.015, '10%', fontsize=8,
          color='#EB8F59', zorder=2)
fig1.scatter(tGAHU2, D17_tGAHU[2], color='#DA6B4C', marker='o',
             zorder=4.5) # 20 % MRm
fig1.text(tGAHU2 + 0.0004, D17_tGAHU[2] - 0.015, '20%', fontsize=8,
          color='#DA6B4C', zorder=4.5)
fig1.scatter(tGAHU3, D17_tGAHU[3], color='#C84840', marker='o',
             zorder=4.5) # 30 % MRm
fig1.text(tGAHU3 + 0.0004, D17_tGAHU[3] - 0.015, '30%', fontsize=8,
          color='#C84840', zorder=4.5)
fig1.scatter(tGAHU4, D17_tGAHU[4], color='#B72433', marker='o',
             zorder=4.5) # 40 % MRm
fig1.text(tGAHU4 + 0.0004, D17_tGAHU[4] - 0.015, '40%', fontsize=8,
          color='#B72433', zorder=4.5)
fig1.scatter(tGAHU5, D17_tGAHU[5], color='#A50026', marker='o',
             zorder=4.5, label='Hebrew University') # 50 % MRm
fig1.text(tGAHU5 + 0.0004, D17_tGAHU[5] - 0.015, '50%', fontsize=8,
          color='#A50026', zorder=4.5)

# x and y labels
fig1.set_xlabel("$\Theta_{\:O_2\:uptake}$")
fig1.set_ylabel("$\Delta'^{17}O$ $O_{2, trop}$ (‰)")

# Legend
fig1.legend(loc='best', edgecolor='#F7F7F7', facecolor='#F7F7F7',
            labelcolor=['#364B9A', '#A50026'])

# Saving figure
plt.tight_layout()
plt.savefig('D17O2vtGA.jpg', dpi=800)

print("The global marine O2 flux of the Mehler reaction assuming it takes up" +
      " 30 percent of gross production is " + str(FO2MRm3) + '.\n')
print("The global marine O2 flux of the Mehler reaction assuming it takes up" +
      " 40 percent of gross production is " + str(FO2MRm4) + '.\n')
print("The global marine O2 flux of the Mehler reaction assuming it takes up" +
      " 50 percent of gross production is " + str(FO2MRm5) + '.')

