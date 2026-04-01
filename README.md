# RIXS EPC Analysis (SrTiO3 and KTaO3)

This repository contains MATLAB code for the analysis of resonant inelastic x-ray scattering (RIXS) spectra to extract mode-resolved electron–phonon coupling (EPC) in SrTiO3 and KTaO3.

The approach is based on multi-phonon modeling of RIXS spectra using externally determined phonon energies.

---

## Overview

The workflow consists of:

1. Visualization of raw RIXS spectra  
2. Multi-phonon fitting of spectra  
3. Extraction of coupling constants (g) and EPC strengths (M)  
4. Convergence testing of intermediate-state cutoffs  
5. Calibration of the elastic line using carbon tape  

---

## Repository Structure

- `main.m`  
  Main analysis script. Performs multi-phonon fitting and extracts EPC parameters.  
  See: `MAIN_README.md`

- `raw_spectra.m`  
  Generates plots of raw RIXS spectra.  
  See: `raw_spectra_README.md`

- `convergence_test.m`  
  Tests convergence of intermediate-state cutoffs.  
  See: `convergence_test_README.md`

- `carbon_tape_fit.m`  
  Determines experimental resolution via elastic peak fitting.  
  See: `carbon_tape_fit_README.md`

- `functions/`  
  Contains material-specific multi-phonon fitting functions.

---

## Data

Experimental HDF5 data are not included.

To run the analysis, data must be placed in:

```
hdf5_data/
```

Access to the data can be requested from the original authors.

---

## Usage

Run main analysis:
```matlab
main
```

Plot raw spectra:
```matlab
raw_spectra
```

Run convergence tests:
```matlab
convergence_test
```

Run carbon tape tests:
```matlab
carbon_tape_fit
```

---

## Acknowledgements

This code is based on earlier implementations developed by Dr. Leonardo Martinelli and Severin Flückiger.

The present version extends these methods for the analysis of RIXS data and the extraction of electron–phonon coupling parameters.