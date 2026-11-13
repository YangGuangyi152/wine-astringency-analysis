# Analysis modules

Run `python run_all.py` for the public workflow. Add `--plots` for figures and `--restricted-dir` for confidential verification.

| Module | Access | Purpose and manuscript use |
|---|---|---|
| `01_ti_processing.py` | Restricted | Gaussian smoothing and the reported TI endpoints; five subquality ratings. Section2.6; Fig.6; TableS3. |
| `02_sensory_statistics.py` | Restricted | Friedman, Kendall W, paired Wilcoxon and within-endpoint Holm correction. Sensory comparisons in Section3.4.1. |
| `03_functional_pca.py` | Public | PCA eigenanalysis, centroids and Hotelling ellipses from aggregate moments. Fig.7A,B. |
| `04_rv_exact_permutation.py` | Public | Direction-encoded RV and all720 label permutations; BH across 156 tests. Fig.7C. |
| `05_hierarchical_clustering.py` | Public | Correlation-distance average-linkage clustering; k=2–5 silhouette screen. Fig.7D. |
| `06_qcmd_processing.py` | Restricted | Recovery summaries from three experiments; saliva stage extraction. Figs.4F,5C and QCM-D methods. |
| `07_qcmd_recovery_statistics.py` | Public | Ordinary one-way ANOVA, Tukey and significance letters from sufficient summaries. Figs.4F,5C. |
| `08_physicochemical_statistics.py` | Public | Bulk endpoint ANOVA/Tukey and stated assumption diagnostics. Fig.2A,C; Section2.7. |
| `09_plot_public_results.py` | Public | Aggregate TI, recovery, fPCA, RV and dendrogram panels. Figs.4–7 numerical panels. |
| `10_plot_restricted_qcmd.py` | Restricted | Representative experiment1 and additional experiments2,3; S2 alignment for additional BSM only. Figs.4B,C,5A and supplementary repeated measurements. |
| `11_fpca_individual_curves.py` | Restricted | Full-curve SVD and individual fPC scores. Fig.7A,B verification. |
| `12_dls_distributions.py` | Public | Original log-grid interpolation and display smoothing. Fig.2E. |
| `13_structure_summary.py` | Public | Mean and sample SD of retained structural indices and total phenolics. Fig.1B–D; Table1. |
| `14_score_metric_correlations.py` | Restricted | Five assessor-level fPC1 Spearman correlations quoted in the text. Section3.4.2. |
| `15_bulk_endpoint_summaries.py` | Public | DLS and electrokinetic summaries for the reported fields. Fig.2C,D; Table2. |
| `16_plot_endpoint_summaries.py` | Public | Structure, bulk, DLS and subquality summary plots. Figs.1,2,6. |
| `verify_results.py` | Public/restricted | Numerical comparison with frozen manuscript targets. Reproducibility check. |
| `common.py; plot_support.py` | Support | Portable paths, CLI, BH adjustment and figure export. Reproducibility support. |

`CODE_MAP.csv` records the inputs and outputs. Only the analyses supporting the manuscript are included. The workflow does not include unrelated exploratory models. TEM images, conceptual diagrams and the Mini Program interface are not produced by these numerical scripts.
