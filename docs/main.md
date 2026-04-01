# main.m: EPC Extraction Workflow

This script performs the full analysis of RIXS spectra to extract electron–phonon coupling (EPC) parameters.

## Overview

For each material (KTO, STO), momentum direction, and scattering angle:

1. Load spectra from HDF5  
2. Preprocess intensity (normalization or background subtraction)  
3. Fit multi-phonon model  
4. Extract coupling constants  
5. Save figures and results  

## Inputs

Defined in the configuration section of the script:

- `phononE`  
  Phonon energies (eV)

- `p0`  
  Initial fit parameters

- `lb`, `ub`  
  Lower and upper parameter bounds

- `fwhm`, `mu`  
  Experimental resolution and pseudo-Voigt mixing

- `anglesDeg`  
  List of scattering angles

## Fit Model

The RIXS spectrum is modeled as a sum of:

- Elastic peak (pseudo-Voigt)  
- Single-phonon contributions  
- Multi-phonon excitations (up to third order)  

The intensity is computed via:

`I(n1,n2,n3,...)`

which evaluates the multi-phonon scattering amplitude including intermediate states.

## Outputs

### Figures

Saved in (generated automatically):

`Figures/`  
`Figures/Elastic_subtracted/`

For each dataset:

- Full fit (including elastic peak)  
- Elastic-subtracted spectrum  
- Individual phonon contributions  

### Tables (MATLAB output)

The script generates tables containing:

- Coupling constants `g1`, `g2`, `g3`  
- EPC strengths `M1`, `M2`, `M3` (meV)  
- 95% confidence intervals  
- Fit residuals  

## Key Functions

- `completeFit_threePhonons`  
  Performs nonlinear least-squares fitting  

- `I(...)`  
  Computes multi-phonon scattering intensity  

- `pseudoVoigtAsymmetric`  
  Models peak shapes  

- `two_theta_to_rlu`  
  Converts scattering angle to momentum transfer  

## Notes

- STO(100) requires background subtraction  
- Elastic width can vary depending on dataset  
- Fit stability depends on initial parameters  

## Typical Runtime

On the order of seconds per spectrum, depending on convergence.

## Purpose

This script extracts mode-resolved electron–phonon coupling from RIXS spectra by fitting the full multi-phonon response.