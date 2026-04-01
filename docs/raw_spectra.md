# raw_spectra.m

This script generates plots of raw RIXS spectra from HDF5 data for all materials and momentum directions.

---

## Overview

For each material and momentum direction:

1. Load RIXS scans from HDF5 files  
2. Extract energy loss and intensity  
3. Apply normalization or background subtraction  
4. Plot spectra for each scattering angle  
5. Save figures  

---

## Data

Input data must be stored in:

```
hdf5_data/
```

Each dataset is expected to contain multiple scans corresponding to different scattering angles.

---

## Outputs

Figures are saved in:

```
Figures/Raw_spectra/<Material>/<Tag>/
```

Each figure corresponds to a single scattering angle.

---

## Notes

- STO(100) spectra require background subtraction  
- Other datasets are normalized by a constant factor  
- Plots are generated with consistent axis limits and formatting
