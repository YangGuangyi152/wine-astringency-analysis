"""Physicochemical endpoint tests and assumption diagnostics."""

from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from common import arguments

ALPHA = 0.05


def tukey(frame):
    t = pairwise_tukeyhsd(frame.value.to_numpy(float), frame.group, alpha=ALPHA)
    rows = []
    for index, (a, b) in enumerate(combinations(t.groupsunique, 2)):
        rows.append(
            dict(
                group_a=a,
                group_b=b,
                difference_b_minus_a=float(t.meandiffs[index]),
                p_adjusted=float(t.pvalues[index]),
                reject=bool(t.reject[index]),
                ci95_low=float(t.confint[index, 0]),
                ci95_high=float(t.confint[index, 1]),
            )
        )
    return pd.DataFrame(rows)


def compact_letters(frame, pairs):
    order = (
        frame.groupby("group").value.mean().sort_values(ascending=False).index.tolist()
    )
    rejected = {
        frozenset([r.group_a, r.group_b]) for r in pairs.itertuples() if r.reject
    }
    cliques = [
        frozenset(c)
        for size in range(1, len(order) + 1)
        for c in combinations(order, size)
        if all((frozenset(pair) not in rejected for pair in combinations(c, 2)))
    ]
    maximal = [c for c in cliques if not any((c < other for other in cliques))]
    maximal.sort(key=lambda c: tuple((order.index(x) for x in order if x in c)))
    labels = {g: "" for g in order}
    for i, clique in enumerate(maximal):
        for g in clique:
            labels[g] += "abcdefghijklmnopqrstuvwxyz"[i]
    for a, b in combinations(order, 2):
        assert bool(set(labels[a]) & set(labels[b])) == (
            frozenset([a, b]) not in rejected
        )
    return labels


def main():
    args = arguments(__doc__)
    raw = pd.read_csv(args.data_dir / "bulk_replicates.csv")
    omnibus = []
    pair_rows = []
    summary_rows = []
    diagnostics = []
    for endpoint, part in raw.groupby("endpoint", sort=False):
        order = [
            g
            for g in ["BSM"] + [f"BSM-F{i}" for i in range(1, 7)]
            if g in set(part.group)
        ]
        part = part[part.group.isin(order)].copy()
        groups = [part.loc[part.group == g, "value"].to_numpy(float) for g in order]
        residual = part.value - part.groupby("group").value.transform("mean")
        shapiro = stats.shapiro(residual)
        levene = stats.levene(*groups, center="median")
        diagnostics.append(
            dict(
                endpoint=endpoint,
                shapiro_W=shapiro.statistic,
                shapiro_p=shapiro.pvalue,
                levene_F=levene.statistic,
                levene_p=levene.pvalue,
            )
        )
        for method in ["ANOVA_Tukey"]:
            f, p = stats.f_oneway(*groups)
            between = len(groups) - 1
            within = sum(map(len, groups)) - len(groups)
            pairs = tukey(part)
            pairs["endpoint"] = endpoint
            pairs["method"] = method
            pair_rows.append(pairs)
            omnibus.append(
                dict(
                    endpoint=endpoint,
                    method=method,
                    F=f,
                    p_value=p,
                    df_between=between,
                    df_within=within,
                )
            )
            s = part.groupby("group").value.agg(["mean", "std", "count"]).reset_index()
            s["endpoint"] = endpoint
            s["method"] = method
            s["letters"] = s.group.map(compact_letters(part, pairs))
            summary_rows.append(s)
    pd.DataFrame(omnibus).to_csv(args.output_dir / "bulk_omnibus.csv", index=False)
    pd.concat(pair_rows, ignore_index=True).to_csv(
        args.output_dir / "bulk_pairwise.csv", index=False
    )
    pd.concat(summary_rows, ignore_index=True).to_csv(
        args.output_dir / "bulk_summary.csv", index=False
    )
    pd.DataFrame(diagnostics).to_csv(
        args.output_dir / "bulk_diagnostics.csv", index=False
    )


if __name__ == "__main__":
    main()
