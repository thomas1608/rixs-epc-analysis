# carbon_tape_fit.m

Fits the elastic peak of a carbon tape reference spectrum to determine the experimental energy resolution.

## Purpose

Used to extract the instrumental FWHM for use in RIXS fitting.

## Workflow

1. Load HDF5 dataset  
2. Select elastic peak window  
3. Fit with:
   - Gaussian  
   - Symmetric pseudo-Voigt  
   - Asymmetric pseudo-Voigt  
4. Compare fits via R²  
5. Save figure  

## Outputs

- Figure:  
  `Figures/Carbon_tape/carbon_tape_fit_comparison.png`

- Printed values:
  - FWHM (Gaussian / pseudo-Voigt)
  - R² values

## Notes

- Asymmetric pseudo-Voigt is typically the most accurate  
- Extracted FWHM is used in the main fitting pipeline  