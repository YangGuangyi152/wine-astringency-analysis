"""Reproduce the five assessor-level fPC1 correlations cited in Section 3.4.2."""

import pandas as pd
from scipy.stats import spearmanr
from common import arguments

METRICS = ["Im_smooth", "TI_AUC_smooth", "duration_above_5_s", "I60", "AUC_60_120"]


def main():
    args = arguments(__doc__, restricted=True)
    scores = pd.read_csv(args.output_dir / "fpca_individual_scores.csv")
    metrics = pd.read_csv(args.output_dir / "participant_metrics.csv").rename(
        columns={"user_id": "assessor_id"}
    )
    merged = scores.merge(
        metrics, on=["assessor_id", "component"], validate="one_to_one"
    )
    rows = []
    for col in METRICS:
        subset = merged[["PC1", col]].dropna()
        rho, p = spearmanr(subset.PC1, subset[col])
        rows.append(
            dict(PC="PC1", metric=col, n=len(subset), spearman_rho=rho, spearman_p=p)
        )
    pd.DataFrame(rows).to_csv(
        args.output_dir / "score_metric_spearman_correlations.csv", index=False
    )


if __name__ == "__main__":
    main()
