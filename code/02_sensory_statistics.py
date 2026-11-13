"""Paired sensory statistics: Friedman, Kendall W, Wilcoxon and Holm."""

import itertools
import numpy as np
import pandas as pd
from scipy import stats
from common import arguments, ORDER

CORE_METRICS = [
    "Im_smooth",
    "Tm_smooth_s",
    "TI_AUC_smooth",
    "AUC_60_120",
    "late_AUC_fraction_60_120",
    "I60",
    "I90",
    "half_life_after_peak_s",
    "duration_above_5_s",
    "dry_level",
    "rough_level",
    "shrink_level",
    "gripping_level",
    "coating_level",
]


def holm_adjust(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    mask = np.isfinite(p)
    vals = p[mask]
    order = np.argsort(vals)
    adjusted = np.empty_like(vals)
    running = 0.0
    m = len(vals)
    for rank, idx in enumerate(order):
        running = max(running, vals[idx] * (m - rank))
        adjusted[idx] = min(running, 1.0)
    out[np.where(mask)[0]] = adjusted
    return out.tolist()


def main():
    args = arguments(__doc__, restricted=True)
    panel = pd.read_csv(args.output_dir / "participant_metrics.csv")
    omnibus = []
    pairs = []
    for metric in CORE_METRICS:
        wide = (
            panel.pivot(index="user_id", columns="component", values=metric)
            .reindex(columns=ORDER)
            .dropna()
        )
        statistic, p = stats.friedmanchisquare(*[wide[c].to_numpy() for c in ORDER])
        omnibus.append(
            dict(
                analysis_set="all_13_panelists",
                metric=metric,
                n_panelists=len(wide),
                friedman_chi2=statistic,
                friedman_p=p,
                kendall_W=statistic / (len(wide) * 5),
            )
        )
        local = []
        for a, c in itertools.combinations(ORDER, 2):
            difference = wide[c] - wide[a]
            if np.allclose(difference, 0):
                w, p = 0.0, 1.0
            else:
                w, p = stats.wilcoxon(
                    wide[c],
                    wide[a],
                    zero_method="wilcox",
                    correction=False,
                    method="auto",
                )
            local.append(
                dict(
                    analysis_set="all_13_panelists",
                    metric=metric,
                    comparison=f"{c}_minus_{a}",
                    mean_first=wide[a].mean(),
                    mean_second=wide[c].mean(),
                    paired_mean_diff_second_minus_first=difference.mean(),
                    wilcoxon_stat=w,
                    wilcoxon_p=p,
                )
            )
        for row, adjusted in zip(local, holm_adjust([r["wilcoxon_p"] for r in local])):
            row["wilcoxon_p_holm_within_metric"] = adjusted
        pairs.extend(local)
    pd.DataFrame(omnibus).to_csv(args.output_dir / "sensory_friedman.csv", index=False)
    pd.DataFrame(pairs).to_csv(args.output_dir / "sensory_pairwise.csv", index=False)


if __name__ == "__main__":
    main()
