"""Recovery ANOVA, Tukey comparisons and letters from public sufficient summaries."""

from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats
from common import arguments

COMPONENTS = []


def letters_from_pairs(summary, pairs, uppercase):
    order = summary.sort_values("mean", ascending=False).index.tolist()
    rejected = {
        frozenset((r.group_a, r.group_b)) for r in pairs.itertuples() if r.reject
    }
    cliques = [
        frozenset(group)
        for size in range(1, 7)
        for group in combinations(order, size)
        if all((frozenset(pair) not in rejected for pair in combinations(group, 2)))
    ]
    maximal = [
        group for group in cliques if not any((group < other for other in cliques))
    ]
    maximal.sort(key=lambda group: tuple((order.index(c) for c in order if c in group)))
    required = {frozenset((c,)) for c in order}
    required.update(
        (
            frozenset(pair)
            for pair in combinations(order, 2)
            if frozenset(pair) not in rejected
        )
    )
    covers = [
        {frozenset((c,)) for c in group}
        | {frozenset(pair) for pair in combinations(group, 2)}
        for group in maximal
    ]
    chosen = None
    for size in range(1, len(maximal) + 1):
        for indices in combinations(range(len(maximal)), size):
            if set().union(*(covers[i] for i in indices)) >= required:
                chosen = indices
                break
        if chosen is not None:
            break
    assert chosen is not None
    alphabet = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if uppercase else "abcdefghijklmnopqrstuvwxyz"
    )
    labels = {c: "" for c in COMPONENTS}
    for index, column in enumerate(chosen):
        for c in maximal[column]:
            labels[c] += alphabet[index]
    for a, b in combinations(COMPONENTS, 2):
        assert bool(set(labels[a]) & set(labels[b])) == (
            frozenset((a, b)) not in rejected
        )
    return labels


def analyse(summary):
    global COMPONENTS
    all_summary = []
    omnibus = []
    all_pairs = []
    for (system, metric, stage), part in summary.groupby(
        ["film_system", "metric", "stage"], sort=False
    ):
        COMPONENTS = sorted(part.component.unique())
        g = part.set_index("component").loc[COMPONENTS]
        means = g["mean"].to_numpy()
        sd = g["std"].to_numpy()
        n = g["count"].to_numpy()
        k = len(g)
        total = int(n.sum())
        grand = np.sum(n * means) / total
        within = np.sum((n - 1) * sd**2)
        between = np.sum(n * (means - grand) ** 2)
        mse = within / (total - k)
        f = (between / (k - 1)) / mse
        omnibus.append(
            dict(
                film_system=system,
                metric=metric,
                stage=stage,
                F=f,
                df_between=k - 1,
                df_within=total - k,
                p_value=stats.f.sf(f, k - 1, total - k),
            )
        )
        pairs = []
        for i, j in combinations(range(k), 2):
            se = np.sqrt(mse / 2 * (1 / n[i] + 1 / n[j]))
            difference = means[j] - means[i]
            p = stats.studentized_range.sf(abs(difference) / se, k, total - k)
            half = stats.studentized_range.ppf(0.95, k, total - k) * se
            pairs.append(
                dict(
                    film_system=system,
                    metric=metric,
                    stage=stage,
                    group_a=COMPONENTS[i],
                    group_b=COMPONENTS[j],
                    mean_difference_b_minus_a=difference,
                    p_adjusted=p,
                    reject=bool(p < 0.05),
                    ci95_low=difference - half,
                    ci95_high=difference + half,
                )
            )
        pairs = pd.DataFrame(pairs)
        g["letters"] = pd.Series(letters_from_pairs(g, pairs, stage == "S3b"))
        all_summary.append(g.reset_index())
        all_pairs.append(pairs)
    return (
        pd.concat(all_summary, ignore_index=True),
        pd.DataFrame(omnibus),
        pd.concat(all_pairs, ignore_index=True),
    )


def main():
    args = arguments(__doc__)
    data = pd.read_csv(args.data_dir / "qcmd_recovery_summary.csv")
    summary, anova, pairs = analyse(data)
    for system in ["BSM", "saliva"]:
        for name, frame in [
            ("recovery_summary", summary),
            ("recovery_anova", anova),
            ("recovery_tukey", pairs),
        ]:
            frame[frame.film_system == system].drop(columns="film_system").to_csv(
                args.output_dir / f"{system}_{name}.csv", index=False
            )


if __name__ == "__main__":
    main()
