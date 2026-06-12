# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains Python code for the chemical reaction network box models and calculations from Valerio 2021 ("Large contribution of light-dependent oxygen uptake to global O2 cycling"), a Rice University thesis on atmospheric oxygen triple-isotope composition.

The models simulate oxygen isotope cycling between stratosphere, troposphere, and biosphere reservoirs to understand factors controlling the Δ'17O signature of atmospheric O2.

## Running the Models

All scripts are standalone Python files. Run directly with Python:
```bash
python Model-Implementations/Valerio-2021-Split-Biosphere-Model/Close-to-Final-Split-Biosphere-Model/Valerio-2021-Close-to-Final-Split-Biosphere-Model.py
```

Model runs are computationally intensive (solving ODEs over 10^6 years with 0.1 year timesteps). Outputs include Excel spreadsheets and matplotlib plots saved to the working directory.

## Dependencies

- numpy
- pandas
- scipy (odeint for ODE integration)
- matplotlib
- seaborn
- labellines (for some plots)

## Code Architecture

### Model Implementations

**Valerio-2021-Split-Biosphere-Model/**: The thesis model that separates terrestrial and marine biospheres into distinct boxes. Key difference from Young et al. 2014 is tracking different respiration/photosynthesis fractionation factors for land vs. ocean.
- `Close-to-Final-Split-Biosphere-Model/`: Production model and figure generation scripts (Figures 6-10)
- `Work-in-Progress-Split-Biosphere-Models/`: Development versions exploring parameter variations

**Young-2014-Python-Implementation/**: Python port of the Young et al. 2014 Fortran model. Combines terrestrial and marine biosphere into one box. Results not used in thesis but useful for validation.

### Model Structure

Each model file follows the same pattern:
1. Initialize moles of isotopologues (Q = 18O, X = 17O notation)
2. Define fractionation factors and rate constants from literature
3. Build system of ~24 coupled ODEs for species in stratosphere/troposphere/biosphere
4. Solve to steady state using scipy.odeint
5. Calculate isotope ratios, delta values, and Δ'17O
6. Compare against target values and export results

### Relevant-Calculations/

Standalone scripts for specific calculations referenced in the thesis (CO2-H2O fractionation factors, marine superoxide fluxes, etc.).

## Key Parameters

- `ft`: Fraction of global primary production from terrestrial biosphere (nominally 0.6)
- `tGA`: Global average oxygen uptake theta (~0.520)
- `alphart`/`alpharm`: Respiration fractionation factors for terrestrial/marine
- `alphaCO2H2O`: CO2-H2O oxygen isotope equilibrium fractionation

## Known Issues

Per the README: the code in this repository may not be the final versions used in the thesis. Some final model calibrations were not pushed before the author lost access to local files.
