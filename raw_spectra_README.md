# raw_spectra.m



Plots raw RIXS spectra from HDF5 files for all materials and momenta.



## Purpose

Visual inspection of spectra before fitting.



## Workflow

1. Load HDF5 scans

2. Extract energy and intensity

3. Normalize or subtract background

4. Plot spectra for each angle

5. Save figures



## Outputs

- Figures in:

&#x20; Figures/Raw\_spectra/<Material>/<Tag>/



## Notes

- STO(100) requires background subtraction
- Other datasets are normalized



