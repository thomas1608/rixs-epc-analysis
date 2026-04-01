\# carbon\_tape\_fit.m



Fits the elastic peak of a carbon tape reference spectrum to determine experimental energy resolution.



\## Purpose

Used to extract the instrumental FWHM for use in RIXS fitting.



\## Workflow

1\. Load HDF5 dataset

2\. Select elastic peak window

3\. Fit with:

&#x20;  - Gaussian

&#x20;  - Symmetric pseudo-Voigt

&#x20;  - Asymmetric pseudo-Voigt

4\. Compare fits via R²

5\. Save figure



\## Outputs

\- Figure: Figures/Carbon\_tape/carbon\_tape\_fit\_comparison.png

\- Printed values:

&#x20; - FWHM (Gaussian / PV)

&#x20; - R² values



\## Notes

\- Asymmetric PV is typically most accurate

\- FWHM used in main fitting pipeline

