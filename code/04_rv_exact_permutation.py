"""Reproduce direction-encoded RV associations and exhaustive six-fraction permutation tests."""

import itertools
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from common import arguments, ORDER, bh_adjust


def rv_coefficient(x, y):
    x = x - x.mean(axis=0, keepdims=True)
    y = y - y.mean(axis=0, keepdims=True)
    xx = x @ x.T
    yy = y @ y.T
    return float(np.sum(xx * yy) / np.sqrt(np.sum(xx * xx) * np.sum(yy * yy)))


def main():
    args = arguments(__doc__)
    manifest = pd.read_csv(args.data_dir / "association_variable_manifest.csv")
    values = pd.read_csv(args.data_dir / "association_component_values.csv")
    vectors = {}
    for col, g in values.groupby("component_column", sort=False):
        x = g.set_index("component").loc[ORDER, "value"].to_numpy(float)
        vectors[col] = ((x - x.mean()) / x.std(ddof=1))[:, None]
    permutations = np.asarray(list(itertools.permutations(range(6))), dtype=int)
    rows = []
    null_rows = []
    for x in manifest.query("role=='experimental'").itertuples():
        for y in manifest.query("role=='sensory'").itertuples():
            xv = vectors[x.component_column]
            yv = vectors[y.component_column]
            r = float(pearsonr(xv[:, 0], yv[:, 0]).statistic)
            rv = rv_coefficient(xv, yv)
            null_r = yv[permutations, 0] @ xv[:, 0] / 5
            null_rv = null_r**2
            exceed = int(np.sum(null_rv >= rv - 1e-12))
            rows.append(
                dict(
                    experimental_order=x.selection_order,
                    experimental_source=x.source,
                    experimental_label=x.display_label,
                    experimental_column=x.component_column,
                    sensory_order=y.selection_order,
                    sensory_label=y.display_label,
                    sensory_column=y.component_column,
                    n_components=6,
                    n_exact_permutations=720,
                    rv_coefficient=rv,
                    rv_exceedance_count=exceed,
                    rv_p_exact_upper_tail=exceed / 720,
                    signed_pearson_r=r,
                    direction_encoded_rv=float(np.sign(r) * rv),
                )
            )
            null_rows.extend(
                dict(
                    experimental_column=x.component_column,
                    sensory_column=y.component_column,
                    permutation_id=i + 1,
                    rv_null=float(value),
                )
                for i, value in enumerate(null_rv)
            )
    result = pd.DataFrame(rows)
    result["rv_q_BH_156"] = bh_adjust(result.rv_p_exact_upper_tail)
    result.to_csv(args.output_dir / "rv_associations.csv", index=False)
    pd.DataFrame(null_rows).to_csv(
        args.output_dir / "exact_null_distributions.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(
        np.asarray(ORDER)[permutations], columns=[f"position_{i+1}" for i in range(6)]
    ).to_csv(args.output_dir / "exact_permutations.csv", index=False)


if __name__ == "__main__":
    main()
