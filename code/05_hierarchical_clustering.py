"""Average-linkage clustering of 12 sensory RV fingerprints."""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
from common import arguments

MODULE_ORDER = ["Peak/load core", "Coarse/pucker", "Persistence bridge", "Late/film"]


def normalize_label(value: str) -> str:
    return " ".join(str(value).split())


def fingerprint_matrix(associations: pd.DataFrame, value_col: str) -> pd.DataFrame:
    sensory_order = (
        associations[["sensory_order", "sensory_label"]]
        .drop_duplicates()
        .sort_values("sensory_order")["sensory_label"]
        .map(normalize_label)
        .tolist()
    )
    table = associations.copy()
    table["sensory_label"] = table["sensory_label"].map(normalize_label)
    table["experimental_label"] = table["experimental_label"].map(normalize_label)
    experimental_order = (
        table[["experimental_order", "experimental_label"]]
        .drop_duplicates()
        .sort_values("experimental_order")["experimental_label"]
        .tolist()
    )
    return table.pivot(
        index="sensory_label", columns="experimental_label", values=value_col
    ).reindex(index=sensory_order, columns=experimental_order)


def correlation_distance(data: np.ndarray) -> np.ndarray:
    distance = squareform(pdist(data, metric="correlation"))
    distance = np.nan_to_num(distance, nan=1.0, posinf=1.0, neginf=1.0)
    np.fill_diagonal(distance, 0)
    return (distance + distance.T) / 2


def assign_modules(
    rv_fingerprint: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    distance = correlation_distance(rv_fingerprint.to_numpy())
    tree = linkage(squareform(distance, checks=False), method="average")
    raw_labels = fcluster(tree, 4, criterion="maxclust")
    raw_members = {
        label: set(rv_fingerprint.index[np.where(raw_labels == label)[0]])
        for label in np.unique(raw_labels)
    }
    module_for_raw: dict[int, str] = {}
    for label, members in raw_members.items():
        if {"Rough", "Pucker"}.issubset(members):
            module_for_raw[label] = "Coarse/pucker"
        elif "Late fraction AUC" in members:
            module_for_raw[label] = "Late/film"
        elif "I90" in members:
            module_for_raw[label] = "Persistence bridge"
        elif "Im" in members:
            module_for_raw[label] = "Peak/load core"
        else:
            raise ValueError(f"Could not identify empirical module: {members}")
    rows = []
    for label, indicator in zip(raw_labels, rv_fingerprint.index):
        module = module_for_raw[int(label)]
        superclass = (
            "Bulk aggregation / peak intensity"
            if module in {"Peak/load core", "Coarse/pucker"}
            else "Interface film / temporal transformation"
        )
        rows.append(
            {
                "sensory_indicator": indicator,
                "raw_cluster_id": int(label),
                "empirical_module": module,
                "mechanism_superclass": superclass,
                "module_order": MODULE_ORDER.index(module) + 1,
            }
        )
    assignments = pd.DataFrame(rows)
    return (assignments, tree, distance)


def linkage_table(tree: np.ndarray, labels: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        tree, columns=["left_node", "right_node", "distance", "n_members"]
    ).assign(
        merge_step=lambda x: np.arange(1, len(x) + 1), leaf_labels=";".join(labels)
    )[
        [
            "merge_step",
            "left_node",
            "right_node",
            "distance",
            "n_members",
            "leaf_labels",
        ]
    ]


def main():
    args = arguments(__doc__)
    associations = pd.read_csv(args.output_dir / "rv_associations.csv")
    fingerprint = fingerprint_matrix(associations, "rv_coefficient")
    assignments, tree, distance = assign_modules(fingerprint)
    assignments.to_csv(args.output_dir / "cluster_assignments.csv", index=False)
    linkage_table(tree, fingerprint.index.tolist()).to_csv(
        args.output_dir / "cluster_linkage.csv", index=False
    )
    fingerprint.to_csv(args.output_dir / "rv_fingerprints.csv")
    pd.DataFrame(distance, index=fingerprint.index, columns=fingerprint.index).to_csv(
        args.output_dir / "clustering_distances.csv"
    )
    rows = []
    for method in ["average"]:
        tree_i = linkage(squareform(distance, checks=False), method=method)
        for k in [2, 3, 4, 5]:
            labels = fcluster(tree_i, k, criterion="maxclust")
            rows.append(
                dict(
                    data_basis="mechanism_fingerprint_rv_coefficient",
                    linkage=method,
                    k=k,
                    silhouette=silhouette_score(distance, labels, metric="precomputed"),
                )
            )
    pd.DataFrame(rows).to_csv(
        args.output_dir / "clustering_parameter_screen.csv", index=False
    )


if __name__ == "__main__":
    main()
