# Data dictionary

CSV files use UTF-8. Blank numeric cells represent unavailable measurements. F1–F6 identify the six tannin fractions. SD uses `ddof=1`. `DATA_SCHEMA.csv` inventories every public CSV field.

## Public data

| File | Scientific content and use |
|---|---|
| `association_component_values.csv` | 150 values:13 experimental and12 sensory variables ×six fractions. Original six-fraction input for RV/clustering (Fig.7C,D). |
| `association_variable_manifest.csv` | Variable roles, labels, selection order and source classification. |
| `fpca_aggregate_moments.npz` | Numeric arrays `time_s` (600), `covariance` (600× 600), `component_centered_mean` (6× 600), `component_n` (six13s), `n_observations` (78). Aggregate moments for fPCA eigenstructure and group centroids. |
| `fpca_group_score_covariance.csv` | Group n and the three distinct entries of each 2×2 score covariance matrix, used for Hotelling confidence ellipses. |
| `sensory_mean_curves.csv` | Six fractions × 600 times, with pointwise mean, SD and count. Intensity0–10, time in seconds (Fig.6). |
| `sensory_component_summary.csv` | Mean and SD of the reported TI endpoints and five subqualities;13 assessors per fraction. |
| `qcmd_recovery_summary.csv` | BSM/saliva ×frequency/dissipation ×S3a/S3b ×fraction: mean, SD and n=3 independent experiments. Units are percentage points (Figs.4F,5C). |
| `qcmd_processed_stage_tables/` | Representative equilibrium and stage-response coordinates used in Figs.4,5. Frequency in Hz, dissipation in ppm, time in minutes. Within-stage point SD is distinguished from between-experiment SD. |
| `structure_replicates.csv` |18 measurements retaining mDP, extension-unit percentage, PD percentage, galloylation percentage and PC/PD. These are LCMS-derived indices, not raw spectra. |
| `structure_summary.csv` | Six fraction means/SD, independently recalculated by script13. |
| `total_phenolics_replicates.csv`; `total_phenolics_summary.csv` | Three observations per fraction and their means/SD; mg GAE/mg (Table1). |
| `bulk_replicates.csv` | A600 absorbance and initial apparent hydrodynamic diameter in nm; BSM and six BSM–tannin mixtures, three measurements/group (Fig.2A,C). |
| `dls_summary_replicates.csv` |21 initial records; sample label, effective diameter, polydispersity, D10, D50 and D90. Diameters in nm (Fig.2C,D). |
| `dls_size_distributions.csv` | Initial/3h distributions for BSM–F1…F6, three measurements; diameter in nm and instrument intensity ordinate (Fig.2E). |
| `zeta_potential_replicates.csv` | Nine valid records for BSM, BSM–F1 and BSM–F2. Retains zeta potential, mobility, conductance and RMS residual, with units in headers (Table2). Missing F3–F6 measurements are not imputed. |

The archived association vectors and the triplicate recovery summaries are documented inputs for their respective analyses. The code preserves the input dataset used for each reported analysis.

## Restricted data

`sensory/records.jsonl` contains 78 unique assessor–fraction records (13× 6). Fields are `assessor_id`, `component`, `time_s`, `intensity`, `ratings`. Assessors are P01–P13. Rating keys map `dry`→drying, `rough`→rough, `shrink`→pucker, `gripping`→grippy, `coating`→mouthcoating; categories 1–5. The reference blank is not assigned a score of 1.

The confidential submission contains source records only. Numerical validation uses the public reference tables and aggregate moments.

`qcmd/BSM/` has 18 trajectories and `qcmd/saliva/` has nine. File prefixes1,2,3 identify independent experiments. Experiment 1 is the representative measurement; experiments2,3 are the additional measurements. Trajectory columns are `time_s`, `delta_f_Hz`, `delta_D_ppm`. Recovery comparisons use the publicly supplied n/mean/sample-SD summaries. The confidential submission does not include additional processed stage or recovery tables.

## Calculation conventions

- TI samples cover 0.2–120 s at nominal 0.2 s intervals. Gaussian smoothing uses σ = 1 s and nearest-edge handling. No artificial t=0 point is inserted. Trapezoidal integration uses stored time coordinates.
- `Im_smooth`: maximum smoothed intensity; `Tm_smooth_s`: first maximum time. `TI_AUC_smooth`: total recorded AUC. `AUC_60_120`: late AUC. `late_AUC_fraction_60_120`: late/total AUC, stored as a fraction. `I60`, `I90`: interpolated residual intensity. `residual_I60_over_Im`: I60/Im. Half-life is the first post-peak time at or below 0.5 Im minus Tm, or missing if not reached. Duration above 5 counts samples strictly above 5 ×median time step.
- fPCA uses the original pointwise-centered, unscaled 600-point curve matrix. The public covariance and confidential SVD reproduce the same eigenstructure. Component signs are fixed using the largest-absolute-loading value. Confidence regions use the 95% Hotelling factor.
- RV uses sample-z scores across six fractions. Univariate RV=r²; displayed direction is sign(Pearson r)×RV. All 720 label permutations, including identity, are enumerated. Upper-tail P uses the original 1e−12 threshold tolerance. BH adjustment is across 156 tests.
- Clustering uses 12 sensory RV profiles across 13 experimental variables, correlation distance and average linkage. Silhouettes screen k=2–5; the reported solution has four groups.
- Sensory comparisons use Friedman/Kendall W and paired two-sided Wilcoxon, with Holm adjustment over 15 fraction pairs for each retained endpoint. The 14 tested endpoints support 210 pairwise comparisons.
- Recovery is 100 ×(S2−S3)/(S2−S1) per experiment. Published recovery comparisons are recalculated from processed per-experiment inputs or sufficient n/mean/sample-SD summaries. Systems, response metrics and rinse stages are separate ANOVA/Tukey families. Different letters indicate adjusted P < 0.05 within a family.

`reference_results` stores frozen numerical targets. Outputs are calculated from inputs before comparison. Numerical checks use rtol=1e−10 and atol=1e−12. Font rendering may vary by platform.
