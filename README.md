# Wine tannin and astringency analyses

Repository: [wine-astringency-analysis](https://github.com/YangGuangyi152/wine-astringency-analysis)  
Analysis version: 1.0.0

This repository contains the analysis code and processed source data supporting the manuscript figures, tables and reported numerical results. Individual sensory records and full QCM-D traces are provided separately for confidential editorial and peer-review verification.

## Installation

Tested with Python 3.12 on macOS ARM64. No GPU is required. Use the pinned requirements, particularly SciPy 1.14.1, to preserve the Wilcoxon automatic handling of ties and zeros used for the reported results.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate `.venv\Scripts\activate`. Installation time depends on the platform and connection. Figure exports use Arial when available; font substitution can affect layout without changing numerical results.

## Public workflow

```bash
python run_all.py
python run_all.py --plots
```

The default workflow reads `data/public`, recalculates the numerical results, writes tables to `results`, and compares them with frozen targets in `reference_results`. `--plots` adds PDF, SVG, PNG and 600-dpi TIFF exports. Expected values include fPC1 81.1%, fPC2 14.9%, a four-cluster silhouette of 0.790, and 156 RV comparisons with 720 exact permutations per comparison. Per-module runtime is saved as `results/run_times.json`.

## Confidential verification

Editors and reviewers who have received the separate restricted data package can run:

```bash
python run_all.py --restricted-dir /path/to/reviewer_only/data --restricted-output-dir /path/to/reviewer_only/results --plots
```

The sensory input contains 78 assessor–fraction records, corresponding to 13 assessors × 6 fractions. Each QCM-D system has three independent recordings, numbered 1, 2, 3. Experiment 1 supplies the representative main-text traces; experiments 2 and 3 supply the additional traces. BSM includes F1–F6; pooled saliva includes F2, F5 and F6.

Restricted outputs are kept outside the public repository. Individual curves, ratings, endpoint tables and fPC scores must not be committed to the public repository.

## Reproduction scope

The public six-fraction vectors reproduce direction-encoded RV coefficients, exhaustive permutation P values, BH correction and hierarchical clustering. Aggregate covariance inputs reproduce the full-curve PCA eigenstructure, centroids and confidence ellipses without publishing individual curves. Recovery ANOVA and Tukey tests are recalculated from sufficient n/mean/sample-SD summaries. Paired sensory rank tests and individual fPC scores use the confidential dataset.

Recovery statistics use the processed per-experiment recovery tables. The archived six-fraction association matrix is the input used for the reported RV analysis; it is retained as its own analysis dataset. Public figure exports include only content supported by public data; individual sensory points and single-experiment QCM-D traces require confidential inputs.

## Documentation

- `CODE_MAP.md` and `CODE_MAP.csv`: module purposes, inputs, outputs and manuscript use.
- `DATA_DICTIONARY.md` and `DATA_SCHEMA.csv`: data definitions and field inventory.
- `VALIDATION.md`: independent numerical verification.
- `SOFTWARE_SUBMISSION_CHECKLIST.md`: software documentation checklist.
- `LICENSE_STATUS.md`: code/data licence status.

The individual sensory records are restricted because consent did not include public release. Full QCM-D traces are restricted because they are part of ongoing laboratory research. Anonymous verification data are supplied confidentially to editors and reviewers. The public repository does not contain these records.

Licensing information is provided in `LICENSE_STATUS.md`.
