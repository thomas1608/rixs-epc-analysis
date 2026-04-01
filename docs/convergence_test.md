\# convergence\_test.m



Tests convergence of intermediate-state cutoffs (m1, m2, m3) in the multi-phonon model.



\---



\## Purpose



Ensure truncation of phonon sums does not affect computed intensities.



\---



\## Workflow



1\. Define EPC parameters (g1,g2,g3)

2\. Sweep m1, m2, m3 independently

3\. Compute total intensity

4\. Detect plateau (relative tolerance)

5\. Save convergence plots



\---



\## Outputs



\- Figures in:



```

Figures/Convergence\_test/

```



Each plot shows:



\- Channel contributions

\- Total intensity

\- Convergence point



\---



\## Notes



\- Convergence depends on g

\- Typical values: m \~ 40–100

