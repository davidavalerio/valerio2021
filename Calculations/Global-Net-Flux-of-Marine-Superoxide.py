#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:08:19 2021

@author: david
"""

import numpy as np

#%% Global marine superoxide flux using Sutherland2020 data

## Parameters from literature

# From cell addition experiments

# Experimental mean net superoxide production in amol / cell h
n1PCMIT9312 = 0.041 # Prochlorococcus MIT9312
n1PCMED4 = 0.005 # Prochlorococcus MED4
n1PCNATL2A = 0.030 # Prochlorococcus NATL2A
n1PCMIT9313 = 0.029 # Prochlorococcus MIT9313
n1SCWH8102VS = 337 # Synechococcus sp. WH8102 in Vineyard Sound media
n1SCWH8102SS = 23.1 # Synechococcus sp. WH8102 in Sargasso Sea media
n1SCWH7803 = 22.8 # Synechococcus sp. WH7803 in Sargasso Sea
n1SARHTCC1062 = 0.15 # Pelagibacter (SAR11) HTCC1062
n1SARHTCC7211 = 0.16 # Pelagibacter (SAR11) HTCC7211
nPA = 3019 # Phaeocystis antarctica
nGC = 6088 # Geminigera cryophila

# From superoxide decay experiments

# Experimental mean gross superoxide production in amol / cell h
gPCMIT9312 = 0.070 # Prochlorococcus MIT9312
gPCMED4 = 0.007 # Prochlorococcus MED4
gPCNATL2A = 0.061 # Prochlorococcus NATL2A
gPCMIT9313 = 0.091 # Prochlorococcus MIT9313
gSC= 786 # Synechococcus sp. WH8102
gSARHTCC1062 = 0.18 # Pelagibacter (SAR11) HTCC1062
gSARHTCC7211 = 0.22 # Pelagibacter (SAR11) HTCC7211

# Experimental mean net superoxide production in amol / cell h
n2PCMIT9312 = 0.019 # Prochlorococcus MIT9312
n2PCMED4 = 0.003 # Prochlorococcus MED4
n2PCNATL2A = 0.023 # Prochlorococcus NATL2A
n2PCMIT9313 = 0.050 # Prochlorococcus MIT9313
n2SC = 103 # Synechococcus sp. WH8102
n2SARHTCC1062 = 0.15 # Pelagibacter (SAR11) HTCC1062
n2SARHTCC7211 = 0.16 # Pelagibacter (SAR11) HTCC7211

# Global mean cell count of relevant microbes in ocean
cPC = 2.9e27 # Prochlorococcus count from Flombaum2013
cSC = 7e26 # Synechococcus count from Flombaum2013
cSAR = 2.4e28 # SAR11 count from Morris2002

## Flux calculations

# For cell addition experiment values

# Global net flux

# Average mean net superoxide production of Prochlorococcus strains
n1PCavg = np.average(np.array([n1PCMIT9312, n1PCMED4, n1PCNATL2A, n1PCMIT9313]))

# Average mean net superoxide production of Synechococcus strains
n1SCavg = np.average(np.array([n1SCWH8102VS, n1SCWH8102SS, n1SCWH7803]))

# Average mean net superoxide production of SAR11 strains
n1SARavg = np.average(np.array([n1SARHTCC1062, n1SARHTCC7211]))

# Global net flux in units of mol / yr
n1FPC = (n1PCavg * cPC + n1SCavg * cSC + n1SARavg * cSAR) * 1e-18 * 8760

print(" Global net flux in units of mol / yr for cell addition" +
      " parameters is " + str(n1FPC) + ".\n")

# For superoxide decay experiment values

# Global gross flux

# Average mean gross superoxide production of Prochlorococcus strains
gPCavg = np.average(np.array([gPCMIT9312, gPCMED4, gPCNATL2A, gPCMIT9313]))

# Average mean gross superoxide production of SAR11 strains
gSARavg = np.average(np.array([gSARHTCC1062, gSARHTCC7211]))

# Global gross flux in units of mol / yr
gFPC = (gPCavg * cPC + gSC * cSC) * 1e-18 * 8760

print(" Global gross flux in units of mol / yr for superoxide decay" +
      " parameters is " + str(gFPC) + ".\n")

# Global net flux

# Average mean net superoxide production of Prochlorococcus strains
n2PCavg = np.average(np.array([n2PCMIT9312, n2PCMED4, n2PCNATL2A, n2PCMIT9313]))

# Average net mean superoxide production of SAR11 strains
n2SARavg = np.average(np.array([n2SARHTCC1062, n2SARHTCC7211]))

# Global net flux in units of mol / yr
n2FPC = (n2PCavg * cPC + n2SC * cSC + n2SARavg * cSAR) * 1e-18 * 8760

print(" Global net flux in units of mol / yr for superoxide decay" +
      " parameters is " + str(n1FPC) + ".")