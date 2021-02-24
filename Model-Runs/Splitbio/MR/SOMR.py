#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:08:19 2021

@author: david
"""

import numpy as np
import pandas as pd
from uncertainties import ufloat
from uncertainties.unumpy import uarray
#from upandas import u_to_separate
import matplotlib.pyplot as plt
import seaborn as sns

#%% Global marine superoxide flux using Sutherland2020 data

## Parameters from literature

# Experimental mean gross superoxide production in amol / cell h
FPCMIT9312 = 0.070 # Prochlorococcus
FPCMED4 = 0.007 # Prochlorococcus
FPCNATL2A = 0.061 # Prochlorococcus
FPCMIT9313 = 0.091 # Prochlorococcus
FSC= 786 # Synechococcus
FSARHTCC1062 = 0.18 # Pelagibacter (SAR11)
FSARHTCC7211 = 0.22 # Pelagibacter (SAR11)

# Global mean cell count of relevant microbes in ocean
cPC = 2.9e27 # Prochlorococcus count from Flombaum2013
cSC = 7e26 # Synechococcus count from Flombaum2013
cSAR = 2.4e28 # SAR11 count from Giovannoni2017

## Flux calculations

# Average mean gross superoxide production of Prochlorococcus ecotypes
FPCavg = np.average(np.array([FPCMIT9312, FPCMED4, FPCNATL2A, FPCMIT9313]))

# Average mean gross superoxide production of SAR11 ecotypes
FSARavg = np.average(np.array([FSARHTCC1062, FSARHTCC7211]))

# Global flux in units of amol / h
gFPC = FPCavg * cPC + FSC * cSC + FSARavg * cSAR

# Global flux in units of mol / yr
gFPC = gFPC * 1e18 * 8760
