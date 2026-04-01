\# RIXS EPC Analysis (SrTiO$\_3$ and KTaO$\_3$)



This repository contains MATLAB code used to analyze resonant inelastic x-ray scattering (RIXS) data in order to extract mode-resolved electron–phonon coupling (EPC) in SrTiO$\_3$ and KTaO$\_3$.



The analysis combines multi-phonon fitting of RIXS spectra with externally determined phonon energies, enabling the extraction of coupling strengths as a function of momentum.



\---



\## Overview



The workflow consists of:



1\. Extraction and visualization of raw RIXS spectra  

2\. Multi-phonon fitting of spectra to obtain EPC parameters  

3\. Convergence testing of intermediate-state cutoffs  

4\. Calibration of the elastic line using carbon tape measurements  



\---



\## File Structure



\- `main.m`  

&#x20; Main analysis script. Performs multi-phonon fitting of RIXS spectra and extracts EPC parameters.



\- `raw\_spectra.m`  

&#x20; Plots raw RIXS spectra from HDF5 files for all measured momenta.



\- `convergence\_test.m`  

&#x20; Tests convergence of intermediate-state cutoffs ($m\_1$, $m\_2$, $m\_3$) used in the multi-phonon model.



\- `carbon\_tape\_fit.m`  

&#x20; Fits the elastic peak using Gaussian and pseudo-Voigt models to determine the experimental resolution.



\- `functions/`

&#x20; - `fitting\_function\_threePhonons\_STO110.m`  

&#x20; - `fitting\_function\_threePhonons\_KTO110.m`  

&#x20; Core fitting routines for STO and KTO datasets.



\---



\## Data



The raw experimental data (HDF5 files) are \*\*not included\*\* in this repository.



To run the analysis, data must be placed in:



hdf5\_data/



Access to the data can be requested from the original authors.



\---



\## Usage



Run the main analysis: main  

Plot raw spectra: raw\_spectra  

Run convergence tests: convergence\_test

---

## Acknowledgements

This code is based on earlier implementations developed by Severin Flückiger and Dr. Leonardo Martinelli.

The present version was extended and adapted for the analysis of RIXS data and extraction of electron–phonon coupling parameters.