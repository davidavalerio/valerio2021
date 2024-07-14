#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 21:54:30 2021

@author: david
"""

#%% Load packages and determine the initial moles of relevant species

import pandas as pd
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns

# Lists for output of for loop
sol = []
fracfluxsol = []

# s values
k41x = np.linspace(0.1, 10)

for k41i in k41x:

    # Terrestrial fraction of global primary production from Field 1998
    ft = .6
        
    # Initial moles of combined isotopologues for relevant species
    Os0 = 1 # O strat – initial moles, any small value works
    O1Ds0 = 1 # O(1D) strat – initial moles, any small value works
    O2s0 = 1 # O2 strat – initial moles, any small value works
    CO2s0 = 4.8e15 # CO2 strat – initial moles, fixes mixing ratio of CO2
    O3s0 = 1 # O3 strat – initial moles, any small value works
    O2t0 = 1 # O2 trop – initial moles, any small value works
    CO2t0 = 4.8e16 # CO2 trop – intial moles (270 ppm), 400ppm = 7.2e16
    O2b0 = 1.83e19 # O2 bio total – from H2O initial moles, not used when H2O is infinite
    O2g0 = 2e17 # O2 geo – moles available for oxidation by O2 trop
    
    # Number of oxygens in each species
    isok1 = pd.Series(np.array([2, 1, 2, 1, 2, 2, 1, 2, 2]), 
                      index = np.array(['Os', 'O1Ds', 'O2s', 'CO2s', 'O3s', 'O2t',
                                        'CO2t','O2b', 'O2g']))
    
    # Multiplier for converting atomic ratios to isotopologue ratios for each species
    isok2 = pd.Series(np.array([1, 1, 2, 1, 3, 2, 1, 2, 2]), 
                      index = np.array(['Os', 'O1Ds', 'O2s', 'CO2s', 'O3s', 'O2t',
                                        'CO2t','O2b', 'O2g']))
    
    # Mole fraction/isotopic ratio used in Young 2014 Fortran Code 
    
    # Q represents 18O, X represents 17O
    rQ = 0.002044928
    rX = 0.000370894
    
    # SMOW
    # =============================================================================
    # rQ = 0.0020052
    # rX = 0.0003799
    # =============================================================================
    
    # Fractional abundance of 18O and 17O from isotope ratios
    def frac(R1, R2):
        rZ = R1 / (1 + R1 + R2)
        return rZ
    
    # Calculate 18O/16O and 17O/16O ratios for each species
    # Z is either 18O (Q) or 17O (X)
    def xxZ(isok2, rZ):
        xxZ = (isok2 * rZ) / (isok2 * rQ + isok2 * rX + 1)
        return xxZ
    
    for species in isok2:
        xxQ = xxZ(isok2, rQ)
        xxX = xxZ(isok2, rX)
        
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
    OQg0 = xxQ[8] * O2g0 # OQ geo - initial moles
    OXg0 = xxX[8] * O2g0 # OX geo - initial moles
    
    # Initial moles of common isotopologues (i.e. no Q or X)
    Os0 = Os0 - Qs0 - Xs0 # O strat - initial moles
    O1Ds0 = O1Ds0 - Q1Ds0 - X1Ds0 # O1D strat - initial moles
    O2s0 = O2s0 - OQs0 - OXs0 # O2 strat - initial moles
    CO2s0 = CO2s0 - COQs0 - COXs0 # CO2 strat - inital moles
    O3s0 = O3s0 - OOQs0 - OOXs0 # O3 strat - initial moles
    O2t0 = O2t0 - OQt0 - OXt0 # O2 strat - initial moles
    CO2t0 = CO2t0 - COQt0 - COXt0 # CO2 strat - initial moles
    O2b0 = O2b0 - OQb0 - OXb0 # O2 bio - initial moles
    O2g0 = O2g0 - OQg0 - OXg0 # O2 geo - initial moles
    
    #%% Constants relevant to calculations
    
    # Fractionation factors
    tequil = 0.528 # nominal TOI equilibration slope
    tevap0 = 0.520 # original evapotranspiration theta from Young2014
    isowater = 1.00525 # tuned photosynthesis source water composition from Young 2014
    alphari = 1 / 1.0182 # isotopic fractionation due to respiration
    #alphar = 1 / isowater * alphari # Calculated base fractionation
    #alphar = 0.9770279689 # from Young 2014 Fortran code 
    #alphart = ft * alphar # original from Young2014
    #alpharm = (1 - ft) * alphar # original from Young2014
    alphaCO2H2O = 1.0413 # from Beck et al. 2005
    
    # Updated fractionation factors
    tresp = 0.5149 # global average oxygen uptake theta consistent with Wostbrock2020
    tevap = 0.519 # evapotranspiration theta from Landais2006
    tphoto = 0.525 # photosynthetic theta from LuzBarkan2011
    isoevap = 1.006825 # isotopic enrichment of water from evapotranspiration from Young2014
    isophoto = 1.0029 # isotopic enrichment of photosynthetic O2 from Eisenstadt 2011
    alphart = 1 / isoevap * alphari  # terrestrial respiration fractionation factor
    alpharm = 1 / isophoto * alphari # marine respiration fractionation factor
    
    # Moles of air in atmosphere from Young 2014
    airs = 1.8e19 # moles of air in the stratosphere
    airt = 1.8e20 # moles of air in the troposphere is about 1/10 that of the stratosphere
    
    # Volume of stratosphere and troposphere in cm^3
    vs = 2.8e25 # volume of stratosphere
    vt = vs / 10 # volume of troposphere is 1/10 that of the stratosphere
    
    # Constants for calculating reaction rates
    secyear = 31556952 # number of seconds in a year
    avo = 6.0221409e23 # Avogadro's number
    
    # Molar masses of oxygen isotopes and carbon from PubChem
    mO = 15.994915 # molar mass of oxygen
    mX = 16.999131 # molar mass of oxygen-17
    mQ = 17.99916  # molar mass of oxygen-18
    mC = 12 # molar mass of carbon
    
    #%% Assign reaction rates  from initial spreadsheet
    
    # Equations used to calculate reaction rates
    
    # Reduced mass equation
    def rm(m1, m2):
        rm = (m1 * m2) / (m1 + m2)
        return rm
    
    # Convert reaction rate from units of cm^3/s to moles/yr
    def tomol(cm3s):
        tomol = cm3s * secyear * (avo / vs)
        return tomol
    
    # Rate constants for transport between boxes
    # 1 - troposphere, 2 - biosphere/hydrosphere, 3 - geosphere, 4 - stratosphere
    # k12 means 1 (trop) -> 2 (bio/hydro) and similarly for others
    k12 = 0.0008 # respiration rate constant yr^-1
    k12t = ft * k12 # respiration rate constant from bio terrestrial yr^-1
    k12m = (1 - ft) * k12 # respiration rate constant from bio marine yr^-1
    k21 = 0.00165 # photosynthesis rate constant yr^-1
    k21t = ft * k21 # photosynthesis rate constant from bio terrestrial yr^-1
    k21m = (1 - ft) * k21 # photosynthesis rate constant from bio marine yr^-1
    k13 = 6e-07 # oxidation rate constant yr^-1
    k31 = 5e-05 # organic burial rate constant yr^-1
    k23 = 1.75e-05 # organic detritus delivery from biosphere to oceans yr^-1
    k23t = ft * k23 # organic detritus delivery from terrestrial biosphere to oceans yr^-1
    k23m = (1 - ft) * k23 # organic detritus delivery from marine biosphere to oceans yr^-1
    k41 = k41i # stat-trop mixing rate constant yr^-1
    k14 = k41 / 10 # trop-strat mixing rate constant yr^-1
    
    # Rate constants for all relevant reactions in mol/yr from initial spreadsheet
    KMIF = 1.065 # MIF for O3 formation
    K1 = 3e-5 # O2 + PHO -> O + O mol/yr
    K2 = K1 # OQ + PHO -> Q + O mol/yr
    K3 = K1 # OX + PHO -> X + O mol/yr
    K4 = 4.6548631e-10 # O2 + O -> O3 mol/yr
    K5 = KMIF * 4.6548631e-10 # O2 + Q -> OOQ mol/yr
    K6 = KMIF * 4.6548631e-10 # O2 + X -> OOX mol/yr
    K7 = KMIF * 4.6548631e-10 # OQ + O -> OOQ mol/yr
    K8 = KMIF * 4.6548631e-10 # OX + O -> OOX mol/yr
    K9 = 1000 # O3 + PHO -> O2 + O mol/yr
    K10 = K9 * (1/3) # OOQ + PHO -> O2 + mol/yr
    K11 = K9 * (2/3) # OOQ + PHO -> OQ + O mol/yr
    K12 = K9 * (1/3) # OOX + PHO -> O2 + X mol/yr
    K13 = K9 * (2/3) # OOX + PHO -> OX + O mol/yr
    K14 = 15768 # O3 + PHO -> O2 + O1D mol/yr
    K15 = K14 * (1/3) # OOQ + PHO -> O2 + Q1D mol/yr
    K16 = K14 * (2/3) # OOQ + PHO -> OQ + O1D mol/yr
    K17 = K14 * (1/3) # OOX + PHO -> O2 + X1D mol/yr
    K18 = K14 * (2/3) # OOX + PHO -> OX + O1D mol/yr
    K19 = 4.6548631e-10 # O3 + O -> O2 + O2 mol/yr
    K20 = K19 # OOQ + O -> O2 + OQ mol/yr
    K21 = K19 # OOX + O -> O2 + OX mol/yr
    K22 = K19 # O3 + Q -> O2 + OQ mol/yr
    K23 = K19 # O3 + X -> O2 + OX mol/yr
    K24 = 8.139058e-05 # O3 + O1D -> O2 + O2 mol/yr
    K25 = K24  # OOQ + O1D -> O2 + OQ mol/yr
    K26 = K24  # OOX + O1D -> O2 + OX mol/yr
    K27 = K24  # O3 + Q1D -> O2 + OQ mol/yr
    K28 = K24  # O3 + X1D -> O2 + OX mol/yr
    K29 = K24  # O3 + O1D -> O2 + O + O mol/yr
    K30 = K24 * (1/2) # OOQ + O1D -> O2 + O + Q mol/yr
    K31 = K24 * (1/2) # OOQ + O1D -> OQ + O + O mol/yr
    K32 = K24 * (1/2) # OOX + O1D -> O2 + O + X mol/yr
    K33 = K24 * (1/2) # OOX + O1D -> OX + O + O mol/yr
    K34 = K24 * (1/2) # O3 + Q1D -> O2 + O + Q mol/yr
    K35 = K24 * (1/2) # O3 + Q1D -> OQ + O + O mol/yr
    K36 = K24 * (1/2) # O3 + X1D -> O2 + O + X mol/yr
    K37 = K24 * (1/2) # O3 + X1D -> OX + O + O mol/yr
    K38 = 1.4e-4 # O2 + O1D -> O2 + O mol/yr
    K39 = K38 # O2 + Q1D -> O2 + Q mol/yr
    K40 = 0 # O2 + Q1D -> OQ + O mol/yr
    K41 = K38 # O2 + X1D -> O2 + X mol/yr
    K42 = 0 # O2 + X1D -> OX + O mol/yr
    K43 = 1.2e-10 # O2 + Q -> OQ + O mol/yr
    K43m = K43 / np.sqrt(rm(2 * mO, mO) / rm(mO + mQ, mO))
    K44 = K43 * (1/2) # OQ + O -> O2 + Q mol/yr
    K45 = K43 # O2 + X -> OX + O mol/yr
    K46 = K43 * (1/2) # OX + O -> O2 + X mol/yr
    K47 = 3.022982e-05 # CO2 + O1D -> CO2 + O mol/yr
    K48 = K47 * (1/2) * (2/3) # CO2 + Q1D -> CO2 + Q mol/yr
    K49 = K47 * (1/2) * (4/3) # CO2 + Q1D -> COQ + O mol/yr
    K50 = K47 * (1/2) * (4/3) # COQ + O1D -> CO2 + Q mol/yr
    K51 = K47 * (1/2) * (2/3) # COQ + O1D -> COQ + O mol/yr
    K52 = K47 * (1/2) * (2/3) # CO2 + X1D -> CO2 + X mol/yr
    K53 = K47 * (1/2) * (4/3) # CO2 + X1D -> COX + O mol/yr
    K54 = K47 * (1/2) * (4/3) # COX + O1D -> CO2 + X mol/yr
    K55 = K47 * (1/2) * (2/3) # COX + O1D -> COX + O mol/yr
    
    #%% Calculate reaction rates
    
    # Calculate reaction rates (mol/yr) from reaction rates (cm3/s) stated in paper
    # =============================================================================
    # K19 = tomol(6.86e-15) # O3 + O -> O2 + O2 1/(yr mol)
    # K20 = K19 * np.sqrt(rm(3 * mO, mO) / rm(2 * mO + mQ, mO)) # OOQ + O -> O2 + OQ 1/(yr mol)
    # K21 = K19 * np.sqrt(rm(3 * mO, mO) / rm(2 * mO + mX, mO)) # OOX + O -> O2 + OX 1/(yr mol)
    # K22 = K19 * np.sqrt(rm(3 * mO, mO) / rm(3 * mO, mQ)) # O3 + Q -> O2 + OQ 1/(yr mol)
    # K23 = K19 * np.sqrt(rm(3 * mO, mO) / rm(3 * mO, mX)) # O3 + X -> O2 + OX 1/(yr mol)
    # K24 = tomol(1.2e-10) # O3 + O1D -> O2 + O2 1/(yr mol)
    # K25 = K24 * np.sqrt(rm(3 * mO, mO) / rm(2 * mO + mQ, mO)) # OOQ + O1D -> O2 + OQ 1/(yr mol)
    # K26 = K24 * np.sqrt(rm(3 * mO, mO) / rm(2 * mO + mX, mO)) # OOX + O1D -> O2 + OX 1/(yr mol)
    # K27 = K24 * np.sqrt(rm(3 * mO, mO) / rm(3 * mO, mQ)) # O3 + Q1D -> O2 + OQ 1/(yr mol)
    # K28 = K24 * np.sqrt(rm(3 * mO, mO) / rm(3 * mO, mX)) # O3 + X1D -> O2 + OX 1/(yr mol)
    # K29 = K24 # O3 + O1D -> O2 + O + O 1/(yr mol)
    # K30 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(2 * mO + mQ, mO)) # OOQ + O1D -> O2 + O + Q 1/(yr mol)
    # K31 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(2 * mO + mQ, mO)) # OOQ + O1D -> OQ + O + O 1/(yr mol)
    # K32 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(2 * mO + mX, mO)) # OOX + O1D -> O2 + O + X 1/(yr mol)
    # K33 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(2 * mO + mX, mO)) # OOX + O1D -> OX + O + O 1/(yr mol)
    # K34 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(3 * mO, mQ)) # O3 + Q1D -> O2 + O + Q 1/(yr mol)
    # K35 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(3 * mO, mQ)) # O3 + Q1D -> OQ + O + O 1/(yr mol)
    # K36 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(3 * mO, mX)) # O3 + X1D -> O2 + O + X 1/(yr mol)
    # K37 = .5 * K24 * np.sqrt(rm(3 *mO, mO) / rm(3 * mO, mX)) # O3 + X1D -> OX + O + O 1/(yr mol)
    # K38 = 1.4e-4 # O2 + O1D -> O2 + O mol/yr
    # K39 = K38 # O2 + Q1D -> O2 + Q mol/yr
    # K40 = 0 # O2 + Q1D -> OQ + O mol/yr
    # K41 = K38 # O2 + X1D -> O2 + X mol/yr
    # K42 = 0 # O2 + X1D -> OX + O mol/yr
    # K43 = tomol(2e-16) * np.sqrt(rm(2 * mO, mO) / rm(mO + mQ, mO)) # O2 + Q -> OQ + O 1/(yr mol)
    # K44 = .5 * tomol(2e-16) * np.sqrt(rm(2 * mO, mO) / rm(2 * mO, mQ)) # OQ + O -> O2 + Q 1/(yr mol)
    # K45 = tomol(2e-16) * np.sqrt(rm(2 * mO, mO) / rm(mO + mX, mO)) # O2 + X -> OX + O 1/(yr mol)
    # K46 = .5 * tomol(2e-16) * np.sqrt(rm(2 * mO, mO) / rm(2 * mO, mX)) # OX + O -> O2 + X 1/(yr mol)
    # K47 = tomol(4.46e-11) # CO2 + O1D -> CO2 + O 1/(yr mol)
    # K48 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + 2 * mO), mQ)) # CO2 + Q1D -> CO2 + Q 1/(yr mol)
    # K49 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + 2 * mO), mQ)) # CO2 + Q1D -> COQ + O 1/(yr mol)
    # K50 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + mO + mQ), mO)) # COQ + O1D -> CO2 + Q 1/(yr mol)
    # K51 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + mO + mQ), mO)) # COQ + O1D -> COQ + O 1/(yr mol)
    # K52 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + 2 * mO), mX)) # CO2 + X1D -> CO2 + X 1/(yr mol)
    # K53 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + 2 * mO), mX)) # CO2 + X1D -> COX + O 1/(yr mol)
    # K54 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + mO + mX), mO)) # COX + O1D -> CO2 + X 1/(yr mol)
    # K55 = .5 * K47 * np.sqrt(rm((mC + 2 * mO), mO) /
    #                          rm((mC + mO + mX), mO)) # COX + O1D -> COX + O 1/(yr mol)
    # =============================================================================
    
    # Hydrosphere rate constants (assuming infinite reservoir)
    kr9 = 0.000930443 # CO2 + H2Q -> COQ + H2O  mol/yr
    kr10 = 0.4370797 # COQ + H2O -> CO2 + H2Q mol/yr
    kr11 = 0.00016875677 # CO2 + H2X -> COX + H2O mol/yr
    kr12 = 0.44544841 # COX + H2O -> CO2 + H2X mol/yr
    
    # Hydrosphere multipliers from original model
    kr9x = kr9 / (alphaCO2H2O * isowater * rQ)
    kr11x = kr11 / (alphaCO2H2O ** tequil * isowater ** tevap0 * rX)
    
    # Splitbio hydrosphere rate constants for base model
    kr10t = ft * kr10 # COQ + H2O -> CO2 + H2Q mol/yr terrestrial
    kr10m = (1 - ft) * kr10 # COQ + H2O -> CO2 + H2Q mol/yr marine
    kr12t = ft * kr12 # COX + H2O -> CO2 + H2X mol/yr terrestrial 
    kr12m = (1 - ft) * kr12 # COX + H2O -> CO2 + H2X mol/yr marine
    
    # Splitbio hydrosphere rate constants for updated model
    kr9t = ft * kr9x * alphaCO2H2O * isoevap * rQ # CO2 + H2Q -> COQ + H2O mol/yr terrestrial
    kr9m = (1 - ft) * kr9x * alphaCO2H2O * isophoto * rQ # CO2 + H2Q -> COQ + H2O mol/yr marine
    kr11t = ft * kr11x  * alphaCO2H2O ** tequil * isoevap ** tevap * rX # CO2 + H2X -> COX + H2O mol/yr terrestrial
    kr11m = (1 - ft) * kr11x * alphaCO2H2O ** tequil * isophoto ** tphoto * rX # COX + H2O -> CO2 + H2X mol/yr marine
    
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
        O2gi = y[24]
        OQgi = y[25]
        OXgi = y[26]
        
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
        
        # X (17O) Stratosphere 
        kXsI = (K3 * OXsi + K12 * OOXsi + K32 * OOXsi * O1Dsi +
                K36 * O3si * X1Dsi + K41 * O2si * X1Dsi + K46 * OXsi * Osi +
                K52 * CO2si * X1Dsi + K54 * COXsi * O1Dsi)
        kXsO = (K6 * O2si + K23 * O3si + K45 * O2si)
        dXs = kXsI - Xsi * kXsO
        
        # Q (18O) Stratosphere
        kQsI = (K2 * OQsi + K10 * OOQsi + K30 * OOQsi * O1Dsi +
                K34 *  O3si * Q1Dsi + K39 * O2si * Q1Dsi + K44 * OQsi * Osi +
                K48 * CO2si * Q1Dsi + K50 * COQsi * O1Dsi)
        kQsO = (K5 * O2si + K22 * O3si + K43 * O2si)
        dQs = kQsI - Qsi * kQsO
        
        # O(1D) Stratosphere
        kO1DsI = (K14 * O3si + K16 * OOQsi + K18 * OOXsi)
        kO1DsO = (K24 * O3si + K25 * OOQsi + K26 * OOXsi + K29 * O3si +
                  K30 * OOQsi + K31 * OOQsi + K32 * OOXsi + K33 * OOXsi +
                  K38 * O2si + K47 * CO2si + K50 * COQsi + K51 * COQsi +
                  K54 * COXsi + K55 * COXsi)
        dO1Ds = kO1DsI - O1Dsi * kO1DsO    
        
        # X(1D) (17O(1D)) Stratosphere
        kX1DsI = (K17 * OOXsi)
        kX1DsO = (K28 * O3si + K36 * O3si + K37 * O3si + K41 * O2si +
                  K42 * O2si + K52 * CO2si + K53 * CO2si)
        dX1Ds = kX1DsI - X1Dsi * kX1DsO
        
        # Q(1D) (18O(1D)) Stratosphere
        kQ1DsI = (K15 * OOQsi)
        kQ1DsO = (K27 * O3si + K34 * O3si + K35 * O3si + K39 * O2si +
                  K40 * O2si + K48 * CO2si + K49 * CO2si)
        dQ1Ds = kQ1DsI - Q1Dsi * kQ1DsO
        
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
        
        # OX (O17O) Stratosphere
        kOXsI = (K13 * OOXsi + K18 * OOXsi + K21 * OOXsi * Osi +
                 K23 * O3si * Xsi + K26 * OOXsi * O1Dsi + K28 * O3si * X1Dsi +
                 K33 * OOXsi * O1Dsi + K37 * O3si * X1Dsi + K42 * O2si * X1Dsi +
                 K45 * O2si * Xsi + k14 * OXti)
        kOXsO = (K3 + K8 * Osi + K46 * Osi + k41)
        dOXs = kOXsI - OXsi * kOXsO
        
        # OQ (O18O) Stratosphere
        kOQsI = (K11 * OOQsi + K16 * OOQsi + K20 * OOQsi * Osi +
                 K22 * O3si * Qsi + K25 * OOQsi * O1Dsi + K27 * O3si * Q1Dsi +
                 K31 * OOQsi * O1Dsi + K35 *  O3si * Q1Dsi + K40 * O2si * Q1Dsi +
                 K43 * O2si * Qsi + k14 * OQti)
        kOQsO = (K2 + K7 * Osi + K44 * Osi + k41)
        dOQs = kOQsI - OQsi * kOQsO
        
        # CO2 Stratosphere
        kCO2sI = (K47 * CO2si * O1Dsi + K48 * CO2si * Q1Dsi + K50 * COQsi * O1Dsi +
                  K52 * CO2si * X1Dsi + K54 * COXsi * O1Dsi + k14 * CO2ti)
        kCO2sO = (K47 * O1Dsi + K48 * Q1Dsi + K49 * Q1Dsi + K52 * X1Dsi +
                  K53 * X1Dsi + k41)
        dCO2s = kCO2sI - CO2si * kCO2sO
        
        # COX (CO17O) Stratosphere
        kCOXsI = (K53 * CO2si * X1Dsi + K55 * COXsi * O1Dsi + k14 * COXti)
        kCOXsO = (K54 * O1Dsi + K55 * O1Dsi + k41)
        dCOXs = kCOXsI - COXsi * kCOXsO
        
        # COQ (CO18O) Stratosphere
        kCOQsI = (K49 * CO2si * Q1Dsi + K51 * COQsi * O1Dsi + k14 * COQti)
        kCOQsO = (K50 * O1Dsi + K51 * O1Dsi + k41)
        dCOQs = kCOQsI - COQsi * kCOQsO
        
        # O3 Stratosphere
        kO3sI = (K4 * O2si * Osi)
        kO3sO = (K9 + K14 + K19 * Osi + K22 * Qsi + K23 * Xsi +
                 K24 * O1Dsi + K27 * Q1Dsi + K28 * X1Dsi + K29 * O1Dsi +
                 K34 * Q1Dsi + K35 * Q1Dsi + K36 * X1Dsi + K37 * X1Dsi)
        dO3s = kO3sI - O3si * kO3sO
        
        # OOX (OO17O) Stratosphere
        kOOXsI = (K6 * O2si * Xsi + K8 * OXsi * Osi)
        kOOXsO = (K12 + K13 + K17 + K18 + K21 * Osi + K26 * O1Dsi +
                  K32 * O1Dsi + K33 * O1Dsi)
        dOOXs = kOOXsI - OOXsi * kOOXsO
        
        # OOQ (OO18O) Stratosphere
        kOOQsI = (K5 * O2si * Qsi + K7 * OQsi * Osi)
        kOOQsO = (K10 + K11 + K15 + K16 + K20 * Osi + K25 * O1Dsi +
                  K30 * O1Dsi + K31 * O1Dsi)
        dOOQs = kOOQsI - OOQsi * kOOQsO
        
        # Troposphere ODEs
        
        # O2 Troposphere
        kO2tI = (k41 * O2si + (k21t + k21m) * O2bi + k31 * O2gi)
        kO2tO = ((k12t + k12m) + k13 + k14)
        dO2t = kO2tI - O2ti * kO2tO
        
        # OX (O17O) Troposphere
        kOXtI = (k41 * OXsi + (k21t + k21m) * OXbi + k31 * OXgi)
        kOXtO = (k12t * (alphart ** tresp) + k12m * (alpharm ** tresp) + k13 + k14)
        dOXt = kOXtI - OXti * kOXtO
        
        # OQ (O18O) Troposphere
        kOQtI = (k41 * OQsi + (k21t + k21m) * OQbi + k31 * OQgi)
        kOQtO = (k12t * alphart + k12m * alpharm + k13 + k14)
        dOQt = kOQtI - OQti * kOQtO
        
        # CO2 Troposphere
        kCO2tI = ((kr10t + kr10m) * COQti + (kr12t + kr12m) * COXti + k41 * CO2si)
        kCO2tO = ((kr9t + kr9m) + (kr11t + kr11m) + k14)
        dCO2t = kCO2tI - CO2ti * kCO2tO
        
        # COX (CO17O) Troposphere
        kCOXtI = ((kr11t + kr11m) * CO2ti + k41 * COXsi) 
        kCOXtO = ((kr12t + kr12m) + k14)
        dCOXt = kCOXtI - COXti * kCOXtO
        
        # COQ (CO18O) Troposphere
        kCOQtI = ((kr9t + kr9m) * CO2ti  + k41 * COQsi)
        kCOQtO = ((kr10t + kr10m) + k14) 
        dCOQt = kCOQtI - COQti * kCOQtO
        
        # Biosphere ODEs (equations not used because we assume an infinite reservoir)
        
       # O2 Biosphere
        # kO2bI = (k12 * O2ti)
        # kO2bO = (k21 + k23)
        # dO2b = kO2bI - O2bi *  kO2bO
        dO2b = 0
        
        # OX (17OO) Biosphere
        # kOXbI = (k12 * ((alphar ** tresp) * OXti)
        # kOXbO = (k21 + k23)
        # dOXb = kOXbI - OXbi * kOXbO
        dOXb = 0
        
        # OQ (18OO) Biosphere
        # kOQbI = (k12 * alphar * O2ti)
        # kOQbO = (k21 + k23)
        # dOQb = kOQbI - OQbi * kOQbO
        dOQb = 0
        
        # Geosphere ODEs
        
        # O2 Geosphere
        kO2gI = (k23t * O2bi + k23m * O2bi + k13 * O2ti)
        kO2gO = k31
        dO2g = kO2gI - O2gi * kO2gO
        
         # OX (O17O) Geosphere
        kOXgI = (k23t * OXbi + k23m * OXbi + k13 * OXti)
        kOXgO = k31
        dOXg = kOXgI - OXgi * kOXgO
        
        # OQ (O18O) Geosphere
        kOQgI = (k23t * OQbi + k23m * OQbi + k13 * OQti)
        kOQgO = k31
        dOQg = kOQgI - OQgi * kOQgO
        
        return np.array([dOs, dXs, dO1Ds, dQs, dX1Ds, dQ1Ds, dO2s, dOXs, dOQs,
                         dCO2s, dCOXs, dCOQs, dO3s, dOOXs, dOOQs, dO2t, dOQt,
                         dOXt, dCO2t, dCOQt, dCOXt,dO2b, dOQb, dOXb, dO2g, dOQg,
                         dOXg])
    
    # Initial conditions for solver
    y0 = np.array([Os0, Xs0, O1Ds0, Qs0, X1Ds0, Q1Ds0, O2s0, OXs0, OQs0, CO2s0,
                   COXs0, COQs0, O3s0, OOXs0, OOQs0, O2t0, OQt0, OXt0, CO2t0,
                   COQt0, COXt0, O2b0, OQb0, OXb0, O2g0,
                   OQg0, OXg0])
    
    
    # Time grid
    t = np.arange(0, 10e5, 0.1)
    
    # Order of species in molar output
    moleso = np.array(['Os', 'Xs', 'O1Ds', 'Qs', 'X1Ds', 'Q1Ds', 'O2s', 'OXs',
                       'OQs', 'CO2s', 'COXs', 'COQs', 'O3s', 'OOXs', 'OOQs',
                       'O2t', 'OQt', 'OXt', 'CO2t', 'COQt', 'COXt', 'O2b',
                       'OQb', 'OXb', 'O2g', 'OQg', 'OXg'])
    
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
    R18_O2b = atomZ(isok2[7], moles.loc['OQb'], moles.loc['O2b'])
    R17_O2b = atomZ(isok2[7], moles.loc['OXb'], moles.loc['O2b'])
    R18_O2g = atomZ(isok2[8], moles.loc['OQg'], moles.loc['O2g'])
    R17_O2g = atomZ(isok2[8], moles.loc['OXg'], moles.loc['O2g'])
    
    # Calculate the d17O and d18O of species.
    d18_Os = deltaZ(R18_Os, rQ)
    d17_Os = deltaZ(R17_Os, rX)
    d18_O1Ds = deltaZ(R18_O1Ds, rQ)
    d17_O1Ds = deltaZ(R17_O1Ds, rX)
    d18_O2s = deltaZ(R18_O2s, rQ)
    d17_O2s = deltaZ(R17_O2s, rX)
    d18_CO2s = deltaZ(R18_CO2s, rQ)
    d17_CO2s = deltaZ(R17_CO2s, rX)
    d18_O3s = deltaZ(R18_O3s, rQ)
    d17_O3s = deltaZ(R17_O3s, rX)
    d18_O2t = deltaZ(R18_O2t, rQ)
    d17_O2t = deltaZ(R17_O2t, rX)
    d18_CO2t = deltaZ(R18_CO2t, rQ)
    d17_CO2t = deltaZ(R17_CO2t, rX)
    d18_O2b = deltaZ(R18_O2b, rQ)
    d17_O2b = deltaZ(R17_O2b, rX)
    d18_O2g = deltaZ(R18_O2g, rQ)
    d17_O2g = deltaZ(R17_O2g, rX)
    
    # Calculate D17 at end of model run
    D17_Os = capD(d17_Os, d18_Os)
    D17_O1Ds = capD(d17_O1Ds, d18_O1Ds)
    D17_O2s = capD(d17_O2s, d18_O2s)
    D17_CO2s = capD(d17_CO2s, d18_CO2s)
    D17_O3s = capD(d17_O3s, d18_O3s)
    D17_O2t = capD(d17_O2t, d18_O2t)
    D17_CO2t = capD(d17_CO2t, d18_CO2t)
    D17_O2b = capD(d17_O2b, d18_O2b)
    D17_O2g = capD(d17_O2g, d18_O2g)
    
    # Order of species in isotope output dataframe
    isotopeso = np.array(['d18_Os', 'd18_O1Ds', 'd18_O2s', 'd18_CO2s', 'd18_O3s',
                          'd18_O2t', 'd18_CO2t', 'd18_O2b', 'd18_O2g', 'd17_Os',
                          'd17_O1Ds', 'd17_O2s', 'd17_CO2s', 'd17_O3s', 'd17_O2t',
                          'd17_CO2t', 'd17_O2b', 'd17_O2g', 'D17_Os', 'D17_O1Ds',
                          'D17_O2s', 'D17_CO2s', 'D17_O3s', 'D17_O2t', 'D17_CO2t',
                          'D17_O2b', 'D17_O2g'])
    
    #%% Mole fraction and flux outputs
    
    # Mole fraction of O2 in troposphere
    def xO2(O2, OX, OQ, air):
        xO2 = (O2 + OX + OQ) / air
        return xO2
    
    # Mole fraction of CO2 in troposphere and stratosphere
    def xCO2(CO2, COX, COQ, air):
        xCO2 = (CO2 + COX + COQ) / air
        return xCO2
    
    # Mole fraction of O3 in stratosphere
    def xO3(O3, OOX, OOQ, air):
        xO3 = (O3 + OOX + OOQ) / air
        return xO3
    
    # D17O XO3 molar flux in per mil moles CO2/yr
    def xJ(xCO2s, capDCO2s, airs, k41):
        xJ = (xCO2s * capDCO2s * airs) / (1 / k41)
        return xJ
    
    # Calculate mole fraction of O2 in the troposphere
    xO2 = xO2(moles.loc['O2t'], moles.loc['OXt'], moles.loc['OQt'], airt)
    
    # Calculate mole fraction of CO2 in the troposphere
    xCO2t = xCO2(moles.loc['CO2t'], moles.loc['COXt'], moles.loc['COQt'], airt)
    
    # Calculate mole fraction of CO2 in the stratosphere
    xCO2s = xCO2(moles.loc['CO2s'], moles.loc['COXs'], moles.loc['COQs'], airs)
    
    # Calculate mole fraction of O3 in the stratosphere
    xO3s = xO3(moles.loc['O3s'], moles.loc['OOXs'], moles.loc['OOQs'], airs)
    
    # Calculate D17O flux
    xJ = xJ(xCO2s, D17_CO2s, airs, k41)
    
    # Order of species in mole fraction and isotope flux output
    fracfluxo = np.array(['xO2', 'xCO2t', 'xCO2s', 'xO3s', 'xJ'])
    
    #%% Output dataframes to spreadsheet
    
    # Create empty dataframes for isotope and fraction/flux outputs
    isotopes = pd.DataFrame(np.full((1, 27), 0, dtype=float), columns = isotopeso)
    fracflux = pd.DataFrame(np.full((1, 5), 0, dtype=float), columns = fracfluxo)
    
    # Assign calculated isotope values to dataframe
    isotopes = pd.DataFrame(np.array([d18_Os, d18_O1Ds, d18_O2s, d18_CO2s,
                                      d18_O3s, d18_O2t, d18_CO2t, d18_O2b,
                                      d18_O2g, d17_Os, d17_O1Ds, d17_O2s,
                                      d17_CO2s, d17_O3s, d17_O2t, d17_CO2t,
                                      d17_O2b, d17_O2g, D17_Os, D17_O1Ds,
                                      D17_O2s, D17_CO2s, D17_O3s, D17_O2t,
                                      D17_CO2t, D17_O2b, D17_O2g]),
                            index = isotopeso)
    
    # Assign calculated mole fraction and flux values to dataframe
    fracflux = pd.DataFrame(np.array([xO2, xCO2t, xCO2s, xO3s, xJ]),
                            index =fracfluxo)
    
     # Add isotopes and fracflux outputs
    sol.append(isotopes)
    fracfluxsol.append(fracflux)

#%% Plots

# D17O2t as function of k41
D17O2t_k41x = []
for i in sol:
    D17O2t = i.loc['D17_O2t'].values
    D17O2t_k41x.append(D17O2t)
D17O2t_k41x = np.hstack(D17O2t_k41x)

# Setting up figure parameters
fig1 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig1 = fig1.add_subplot(1, 1, 1)
fig1.set(xlim = (0, 10), ylim = (-0.5, -0.35))

# Plotting D17O2t as function of k41
fig1.plot(k41x, D17O2t_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, -1, -0.459, color='red', linestyles='dotted',
           label = 'k41 = 1 $yr^{-1}$, $^{17}\Delta$ $O_{2, trop}$ = -0.459‰')
plt.hlines(-0.459, 0, 1, color='red', linestyles='dotted')

# Legend and title
fig1.set_xlabel("k41 ($yr^{-1}$)")
fig1.set_ylabel("$^{17}\Delta$ $O_{2, trop}$ (‰)")
fig1.set_title('$^{17}\Delta$ $O_{2, trop}$ vs. k41')
fig1.legend(loc = 'best', fontsize = 'small')
plt.tight_layout()

# Saving plot
plt.savefig('Young2014D17O2tk41.png', dpi = 800)

# d'18O O3 as function of k41
d18O3s_k41x = []
for i in sol:
    d18O3s = i.loc['d18_O3s'].values
    d18O3s_k41x.append(d18O3s)  
    
d18O3s_k41x = np.hstack(d18O3s_k41x)

# Setting up figure parameters
fig2 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig2 = fig2.add_subplot(1, 1, 1)
fig2.set(xlim = (0, 10), ylim = (80, 100))

# Plotting d18O3s as function of k41
fig2.plot(k41x, d18O3s_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, 0, 86, color='red', linestyles='dotted',
           label = "k41 = 1 $yr^{-1}$, $\delta'^{18}$O O$_{3, strat}$ = 86‰")
plt.hlines(86, 0, 1, color='red', linestyles='dotted')

# Plotting target value from Young2014
plt.hlines(94, 0, 10, color='green', linestyles='dashed',
           label = "Target $\delta'^{18}$O O$_{3, strat}$ = 94 ± 13‰")

# Legend and title
fig2.set_xlabel("k41 ($yr^{-1}$)")
fig2.set_ylabel("$\delta'^{18}$O O$_{3, strat}$ (‰)")
fig2.set_title("$\delta'^{18}$O O$_{3, strat}$ vs. k41")
fig2.legend(loc = 'best', fontsize = 'small')

# Saving plot
plt.savefig('Young2014d18O3sk41.png', dpi = 800)
plt.tight_layout()

# d'17O O3 as function of K41
d17O3s_k41x = []
for i in sol:
    d17O3s = i.loc['d17_O3s'].values
    d17O3s_k41x.append(d17O3s)
    
d17O3s_k41x = np.hstack(d17O3s_k41x)

# Setting up figure parameters
fig3 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig3 = fig3.add_subplot(1, 1, 1)
fig3.set(xlim = (0, 10), ylim = (70, 80))

# Plotting d'17O O3 as function of k41
fig3.plot(k41x, d17O3s_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, 0, 74.671, color='red', linestyles='dotted',
           label = "k41 = 1 $yr^{-1}$, $\delta'^{17}$O O$_{3, strat}$ = 74.671‰")
plt.hlines(74.671, 0, 1, color='red', linestyles='dotted')

# Plotting target value from Young2014
plt.hlines(77, 0, 10, color='green', linestyles='dashed',
           label = "Target $\delta'^{17}$O O$_{3, strat}$ = 77 ± 10‰")

# Legend and title
fig3.set_xlabel("k41 ($yr^{-1}$)")
fig3.set_ylabel("$\delta'^{17}$O O$_{3, strat}$ (‰)")
fig3.set_title("$\delta'^{17}$O O$_{3, strat}$ vs. k41")
fig3.legend(loc = 'best', fontsize = 'small')
plt.tight_layout()

# Saving plot
plt.savefig('Young2014d17O3sk41.png', dpi = 800)

# D17O3s as function of k41
D17O3s_k41x = []
for i in sol:
    D17O3s = i.loc['D17_O3s'].values
    D17O3s_k41x.append(D17O3s)
    
D17O3s_k41x = np.hstack(D17O3s_k41x)

# Setting up figure parameters
fig4 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig4 = fig4.add_subplot(1, 1, 1)
fig4.set(xlim = (0, 10), ylim = (25, 35))

# Plotting D17O3s as a function of k41
fig4.plot(k41x, D17O3s_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, -10, 29.262, color='red', linestyles='dotted',
           label = "k41 = 1 $yr^{-1}$, $\Delta'^{17}$O O$_{3, strat}$ = 29.262‰")
plt.hlines(29.262, 0, 1, color='red', linestyles='dotted')

# Plotting target value from Young2014
plt.hlines(29, 0, 10, color='green', linestyles='dashed',
           label = "Target $\Delta'^{17}$ O$_{3, strat}$ = 29 ± 4‰")

# Legend and title
fig4.set_xlabel("k41 ($yr^{-1}$)")
fig4.set_ylabel("$^{17}\Delta$ O$_{3, strat}$ (‰)")
fig4.set_title("$^{17}\Delta$ O$_{3, strat}$ vs. k41")
fig4.legend(loc = 'best', fontsize = 'small')
plt.tight_layout()

# Saving plot
plt.savefig('Young2014D17O3sk41.png', dpi = 800)

# d18O2t as function of k41
d18O2t_k41x = []
for i in sol:
    d18O2t = i.loc['d18_O2t'].values
    d18O2t_k41x.append(d18O2t)  
    
d18O2t_k41x = np.hstack(d18O2t_k41x)

# Setting up figure parameters
fig5 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig5 = fig5.add_subplot(1, 1, 1)
fig5.set(xlim = (0, 10), ylim = (23, 23.3))

# Plotting d18O3s as function of k41
fig5.plot(k41x, d18O2t_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, 0, 23.033, color='red', linestyles='dotted',
           label = "k41 = 1 $yr^{-1}$, $\delta'^{18}$O O$_{2, trop}$ = 23.033‰")
plt.hlines(23.033, 0, 1, color='red', linestyles='dotted')

# Plotting target value from Young2014
plt.hlines(23.2, 0, 10, color='green', linestyles='dashed',
           label = "Target $\delta'^{18}$O O$_{2, trop}$ = 23.2‰")

# Legend and title
fig5.set_xlabel("k41 ($yr^{-1}$)")
fig5.set_ylabel("$\delta'^{18}$O O$_{2, trop}$ (‰)")
fig5.set_title("$\delta'^{18}$O O$_{2, trop}$ vs. k41")
fig5.legend(loc = 'best', fontsize = 'small')
plt.tight_layout()

# Saving plot
plt.savefig('Young2014d18O2tk41.png', dpi = 800)

# d18CO2t as function of k41
d18CO2t_k41x = []
for i in sol:
    d18CO2t = i.loc['d18_CO2t'].values
    d18CO2t_k41x.append(d18CO2t)  
    
d18CO2t_k41x = np.hstack(d18CO2t_k41x)

# Setting up figure parameters
fig6 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig6 = fig6.add_subplot(1, 1, 1)
fig6.set(xlim = (0, 10), ylim = (40, 41.25))

# Plotting d18CO2s as function of k41
fig6.plot(k41x, d18CO2t_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, 0, 40.856, color='red', linestyles='dotted',
           label = "k41 = 1 $yr^{-1}$, $\delta'^{18}$O CO$_{2, trop}$ = 40.856‰")
plt.hlines(40.856, 0, 1, color='red', linestyles='dotted')

# Plotting target value from Young2014
plt.hlines(40.2, 0, 10, color='green', linestyles='dashed',
           label = "Target $\delta'^{18}$O CO$_{2, trop}$ = 40.2‰")

# Legend and title
fig6.set_xlabel("k41 ($yr^{-1}$)")
fig6.set_ylabel("$\delta'^{18}$O CO$_{2, trop}$ (‰)")
fig6.set_title("$\delta'^{18}$O CO$_{2, trop}$ vs. k41")
fig6.legend(loc = 'best', fontsize = 'small')
plt.tight_layout()

# Saving plot
plt.savefig('Young2014d18CO2tk41.png', dpi = 800)

# D17CO2s as function of k41
D17CO2s_k41x = []
for i in sol:
    D17CO2s = i.loc['D17_CO2s'].values
    D17CO2s_k41x.append(D17CO2s)
    
D17CO2s_k41x = np.hstack(D17CO2s_k41x)

# Setting up figure parameters
fig7 = plt.figure(figsize = (5, 5))
with sns.axes_style("whitegrid"):
    fig7 = fig7.add_subplot(1, 1, 1)
fig7.set(xlim = (0, 10), ylim = (0.5, 3))

# Plotting D17O3s as a function of k41
fig7.plot(k41x, D17CO2s_k41x, color='blue')

# Plotting k41 from Young2014
plt.vlines(1, -1, 2.372, color='red', linestyles='dotted',
           label = 'k41 = 1 $yr^{-1}$, $^{17}\Delta$ $CO_{2, strat}$ = 2.372‰')
plt.hlines(2.372, 0, 1, color='red', linestyles='dotted')

# Plotting target value from Young2014
plt.hlines(1.65, 0, 10, color='green', linestyles='dashed',
           label = "Target $^{17}\Delta$ $CO_{2, strat}$ = 1.65‰")

# Legend and title
fig7.set_xlabel("k41 ($yr^{-1}$)")
fig7.set_ylabel("$^{17}\Delta$ $CO_{2, strat}$ (‰)")
fig7.set_title('$^{17}\Delta$ $CO_{2, strat}$ vs. k41')
fig7.legend(loc = 'best', fontsize = 'small')
plt.tight_layout()

# Saving plot
plt.savefig('Young2014D17CO2sk41.png', dpi = 800)