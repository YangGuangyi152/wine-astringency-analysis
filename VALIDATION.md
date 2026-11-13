# Numerical validation

A copy of the release source and inputs was executed in a separate directory. Individual sensory and QCM-D inputs were stored outside the public repository. No calculation input was read from an original project path.

- Public workflow: **18 numerical checks passed**.
- Full confidential workflow: **33 numerical checks passed**.
- The retained 14 sensory endpoints yield 210 paired comparisons; their Holm-adjusted P values reproduce the archived values.
- fPC1/fPC2 variance proportions reproduce 81.1% and 14.9%.
- All 156 RV coefficients, exact permutation P values and BH q values reproduce the reported analysis; each test enumerates 720 permutations.
- The four module memberships and linkage reproduce the archived solution; silhouette rounds to 0.790.
- BSM and saliva recovery means, sample SD, ANOVA, Tukey and significance letters reproduce the documented processed recovery inputs.
- The five assessor-level fPC1 correlations quoted in the text reproduce the archived coefficients and P values.
- Fourteen figure sets passed panel alignment, text-size and collision checks. All effective text sizes are at least 6 pt.

Floating-point comparisons use rtol=1e−10 and atol=1e−12. Fonts and image rasterization may vary by platform; numerical agreement does not require identical output-file bytes. Frozen reference tables are read for validation only, after calculation.
