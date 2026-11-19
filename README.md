# Wine tannin and astringency analyses

Analysis code and processed data for tannin structure, bulk aggregation, QCM-D film responses and time-intensity (TI) sensory analysis.

## Installation and execution

Tested with Python 3.12 on macOS ARM64 and the package versions in `requirements.txt`. No GPU is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
```

On Windows, activate the environment with `.venv\Scripts\activate`.

The default workflow reads `data/public` and writes calculated tables to `results`. To export figures as PDF, SVG, PNG and 600-dpi TIFF:

```bash
python run_all.py --plots
```

An alternative output folder can be specified with `--output-dir /path/to/results`. Figure exports use Arial when available; another font may be substituted. The scripts produce analysis tables and figures without runtime logs or verification reports. Generated outputs are not included in the source package.

## Analysis methods and expected results

- **Functional PCA:** the pointwise-centred, unscaled 600-point TI curves are analysed using their aggregate covariance matrix. The first two components explain 81.1% and 14.9% of the variance. Loading signs are fixed by making the largest absolute loading positive. Group covariance matrices and the 95% Hotelling factor are used for centroid confidence ellipses.
- **RV associations:** 12 experimental variables are compared with 12 sensory indicators across F1–F6, giving 144 comparisons. Sample z scores are used. For univariate inputs, RV = r²; the heatmap displays sign(Pearson r) × RV. All 720 label permutations, including the identity permutation, are enumerated. The upper-tail P value is the proportion with RV at least as large as the observed value, with a numerical tolerance of 1e−12. Benjamini–Hochberg adjustment is applied jointly to all 144 P values; adjusted values are stored in `rv_q_BH_144`.
- **Hierarchical clustering:** each sensory indicator is represented by its 12 unsigned RV coefficients. Correlation distance and average linkage are used, and k = 2–5 is evaluated by silhouette score. The highest score is obtained at k = 4 (0.784549, reported as 0.785). The groups are Peak/load core, Coarse/pucker, Persistence bridge and Late/film. The heatmap and dendrogram use the same sensory-indicator order.
- **QCM-D recovery:** mean, sample SD and n are supplied for three independent experiments. One-way ANOVA and Tukey comparisons are recalculated from these sufficient summaries, separately for each film system, response metric and rinsing stage. Lowercase and uppercase letters denote comparisons within S3a and S3b, respectively; the two stages are not compared with each other.
- **Bulk endpoints:** one-way ANOVA, Tukey comparisons, residual Shapiro–Wilk tests and median-centred Levene tests are calculated from the replicate data.
- **Structural and electrophoretic summaries:** means, sample SDs and counts are calculated from the supplied replicate tables. Missing electrophoretic measurements are not imputed.

## Scripts

Numbered scripts can be run individually with `--data-dir` and `--output-dir`. Scripts that use an earlier script's output must be run in the order below or through `run_all.py`.

| Script in `code/` | Purpose |
|---|---|
| `03_functional_pca.py` | PCA eigenstructure, fraction centroids and confidence ellipses from aggregate moments. |
| `04_rv_exact_permutation.py` | The 144 RV coefficients, exact permutation tests and BH adjustment. |
| `05_hierarchical_clustering.py` | RV fingerprints, distance matrix, linkage, group assignments and silhouette screening. |
| `07_qcmd_recovery_statistics.py` | Recovery ANOVA, Tukey comparisons and significance letters. |
| `08_physicochemical_statistics.py` | Bulk endpoint comparisons and assumption diagnostics. |
| `12_dls_distributions.py` | Log-grid interpolation and smoothing of DLS display curves. |
| `13_structure_summary.py` | Structural indices and total phenolics summaries. |
| `15_bulk_endpoint_summaries.py` | DLS and electrophoretic endpoint summaries. |
| `09_plot_public_results.py` | TI mean curves, recovery, PCA, RV and clustering plots. |
| `16_plot_endpoint_summaries.py` | Structural, bulk, DLS and subquality plots. |
| `01_ti_processing.py` | Individual TI endpoints, ratings and group summaries; confidential input required. |
| `02_sensory_statistics.py` | Friedman/Kendall W and paired Wilcoxon–Holm tests; requires script 01. |
| `06_qcmd_processing.py` | Summaries from optional per-experiment recovery tables; saliva stage means when stage windows are supplied. |
| `11_fpca_individual_curves.py` | PCA from individual curves; requires script 01. |
| `14_score_metric_correlations.py` | Assessor-level fPC1 correlations; requires scripts 01 and 11. |
| `10_plot_restricted_qcmd.py` | Representative and additional QCM-D traces; confidential input required. |
| `common.py`; `plot_support.py` | Shared input handling, BH adjustment and figure exports. |

## Public data

CSV files use UTF-8. F1–F6 identify tannin fractions. Blank numeric entries denote unavailable measurements. SD is calculated with ddof = 1.

| File in `data/public/` | Contents |
|---|---|
| `association_variable_manifest.csv` | Definitions of 24 variables: 12 experimental and 12 sensory. Fields give order, source category, role, display label and variable identifier. |
| `association_component_values.csv` | Six fraction values for each variable: 144 input rows, with fields `component_column`, `component` and `value`. These are the inputs for RV analysis and clustering. |
| `fpca_aggregate_moments.npz` | Numeric arrays: `time_s` (600), `covariance` (600 × 600), `component_centered_mean` (6 × 600), `component_n` (six counts of 13) and `n_observations` (78). |
| `fpca_group_score_covariance.csv` | Fraction counts and the three distinct entries of each 2 × 2 fPC score covariance matrix. |
| `sensory_mean_curves.csv` | Mean, SD and count at each of 600 time points for each fraction. Time in seconds; intensity on a 0–10 scale. |
| `sensory_component_summary.csv` | Fraction means and SDs of TI endpoints and five subquality ratings, with 13 assessors per fraction. |
| `qcmd_recovery_summary.csv` | Frequency/dissipation recovery mean, SD and count for each fraction, system and rinsing stage. Recovery in percent; SD in percentage points. |
| `qcmd_processed_stage_tables/` | Representative stage-response coordinates and equilibrium summaries for BSM and pooled saliva. Frequency in Hz, dissipation in ppm and time in minutes. Within-stage SD describes sampled points, not variation across independent experiments. |
| `structure_replicates.csv`; `structure_summary.csv` | Replicate LC-MS-derived mDP, prodelphinidin percentage, galloylation percentage and procyanidin/prodelphinidin ratio, and their means/SDs. |
| `total_phenolics_replicates.csv`; `total_phenolics_summary.csv` | Total phenolics replicates and means/SDs in mg GAE/mg. |
| `bulk_replicates.csv` | A600 and initial apparent hydrodynamic diameter for BSM and six BSM–tannin mixtures, with three observations per group. Diameter in nm. |
| `dls_summary_replicates.csv` | Initial effective diameter, polydispersity, D10, D50 and D90 for 21 records. Diameters in nm. |
| `dls_size_distributions.csv` | Initial and 3 h DLS distributions; diameter in nm and instrument intensity ordinate. |
| `zeta_potential_replicates.csv` | Nine valid records for BSM, BSM–F1 and BSM–F2, with units in column headers. |

The experimental RV variables comprise mDP, PD%, A600, effective diameter, D90, S2 frequency/dissipation responses and their ratio, and frequency/dissipation recovery after S3a and S3b. The 12 sensory indicators are Im, total AUC, late AUC fraction, I60, I90, post-peak half-life, duration above 5, rough, pucker, drying, grippy and mouthcoating. The association table supplies the six-fraction vectors used for this analysis.

## Confidential source records

Individual sensory records are not public because participant consent did not include public release. Full QCM-D traces are restricted because they form part of ongoing laboratory research. Anonymised source records are supplied separately for confidential editorial and peer review. They are not included in this repository.

The public data reproduce the RV/permutation/BH analysis, clustering, aggregate PCA results, recovery ANOVA/Tukey comparisons and the supplied physicochemical summaries. Individual paired sensory tests, individual fPC scores and raw QCM-D trace plots require the confidential inputs.

Editors and reviewers with that separate package can run:

```bash
python run_all.py --restricted-dir /path/to/reviewer_only/data --restricted-output-dir /path/to/private_results --plots
```

The confidential input layout is `sensory/records.jsonl` and `qcmd/{BSM,saliva}/`. The sensory file contains 78 records (13 assessors × 6 fractions), with fields `assessor_id`, `component`, `time_s`, `intensity` and `ratings`. Assessor codes are P01–P13. Rating keys `dry`, `rough`, `shrink`, `gripping` and `coating` denote drying, rough, pucker, grippy and mouthcoating, respectively, on a five-point scale.

QCM-D trajectories are named `{experiment}_{fraction}.csv`, with columns `time_s`, `delta_f_Hz` and `delta_D_ppm`. Experiments 1–3 are independent; experiment 1 is the representative record. There are 18 BSM traces (F1–F6) and nine pooled-saliva traces (F2, F5 and F6). BSM experiments 2 and 3 are shifted by a constant per trace to align S2 onset for display; this does not change the recorded signal values. Recovery tests use the public recovery summaries. Optional `recovery.csv` and saliva `stage_windows.csv` inputs can be processed by script 06 if supplied separately.

TI processing uses a Gaussian kernel with σ = 1 s and nearest-edge handling. Stored sampling times are 0.2–120 s at nominal 0.2 s intervals; no t = 0 point is inserted. AUC is calculated by trapezoidal integration. Half-life is the first post-peak time at or below half the maximum, minus the peak time; it is left missing if the threshold is not reached. Sensory pairwise tests use SciPy 1.14.1 Wilcoxon `method="auto"` with Holm adjustment across 15 fraction pairs within each endpoint. Confidential outputs must remain outside the public repository.
