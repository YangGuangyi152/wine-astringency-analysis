"""Draw the numerical figure panels reproducible from publicly shared data."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.cluster.hierarchy import dendrogram
from common import arguments, ORDER
from plot_support import style, save, COLORS


def main():
    args = arguments(__doc__)
    style()
    folder = args.output_dir / "figures"
    folder.mkdir(exist_ok=True)
    curves = pd.read_csv(args.data_dir / "sensory_mean_curves.csv")
    fig, ax = plt.subplots(figsize=(7.1, 4.15))
    for c, g in curves.groupby("component", sort=True):
        ax.plot(g.time_s, g["mean"], color=COLORS[c], label=c)
        ax.fill_between(
            g.time_s,
            g["mean"] - g["std"],
            g["mean"] + g["std"],
            color=COLORS[c],
            alpha=0.12,
        )
    ax.set(
        xlabel="Time (s)", ylabel="Astringency intensity", xlim=(0, 120), ylim=(0, 10)
    )
    ax.legend(ncol=3, frameon=False)
    ax.tick_params(axis="x", pad=10)
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.22, top=0.96)
    save(fig, folder / "sensory_mean_TI")
    for system in ["BSM", "saliva"]:
        data = pd.read_csv(args.output_dir / f"{system}_recovery_summary.csv")
        fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.65))
        fig.subplots_adjust(left=0.12, right=0.985, bottom=0.14, top=0.90, wspace=0.34)
        components = sorted(data.component.unique())
        x = np.arange(len(components))
        width = 0.34
        for ax, metric in zip(axes, ["delta_f", "delta_d"]):
            ax.set_title(
                "Frequency recovery" if metric == "delta_f" else "Dissipation recovery"
            )
            ax.axhline(100, color="#BDBDBD", lw=0.85, ls=(0, (3, 3)), zorder=0)
            for offset, stage, alpha, hatch in [
                (-width / 2, "S3a", 0.35, "////////" if system == "saliva" else "////"),
                (width / 2, "S3b", 0.92, ""),
            ]:
                rows = (
                    data[(data.metric == metric) & (data.stage == stage)]
                    .set_index("component")
                    .loc[components]
                )
                ax.bar(
                    x + offset,
                    rows["mean"],
                    width,
                    color=[
                        plt.matplotlib.colors.to_rgba(COLORS[c], alpha)
                        for c in components
                    ],
                    edgecolor="#2A2A2A",
                    linewidth=0.65,
                    hatch=hatch,
                    yerr=rows["std"],
                    capsize=3,
                    error_kw={"elinewidth": 1.0, "capthick": 1.0},
                    rasterized=True,
                )
                for position, row in zip(x + offset, rows.itertuples()):
                    y = row.mean + row.std + 2.5
                    if y < 100 and y + 8 > 100:
                        y = 104
                    ax.text(
                        position, y, row.letters, ha="center", va="bottom", fontsize=10
                    )
            ax.set_xticks(x, components)
            ax.set_ylabel("Recovery (%)", fontweight="bold")
            ax.set_ylim(0, 120)
            ax.set_yticks(np.arange(0, 101, 20))
            ax.set_box_aspect(1)
            (
                ax.set_xlim(-0.48, len(components) - 0.52)
                if system == "saliva"
                else ax.set_xlim(-0.65, 5.65)
            )
        save(fig, folder / f"{system}_recovery")
    ell = pd.read_csv(args.output_dir / "fpca_centroids_and_ellipses.csv")
    loads = pd.read_csv(args.output_dir / "fpca_loadings.csv")
    variance = pd.read_csv(args.output_dir / "fpca_variance.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.17, top=0.86, wspace=0.3)
    for row in ell.itertuples():
        axes[0].add_patch(
            Ellipse(
                (row.PC1_mean, row.PC2_mean),
                row.ellipse_width_95,
                row.ellipse_height_95,
                angle=row.ellipse_angle_deg,
                color=COLORS[row.component],
                alpha=0.18,
            )
        )
        axes[0].plot(
            row.PC1_mean,
            row.PC2_mean,
            "o",
            color=COLORS[row.component],
            label=row.component,
        )
    axes[0].legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=6,
        frameon=False,
        fontsize=8,
        handletextpad=0.25,
        columnspacing=0.8,
    )
    axes[0].autoscale_view()
    axes[0].set_xlabel(f"fPC1 ({variance.explained_variance_ratio.iloc[0]*100:.1f}%)")
    axes[0].set_ylabel(f"fPC2 ({variance.explained_variance_ratio.iloc[1]*100:.1f}%)")
    for pc in ["PC1", "PC2"]:
        axes[1].plot(loads.time_s, loads[f"{pc}_loading"], label=pc)
    axes[1].set(xlabel="Time (s)", ylabel="Loading")
    axes[1].legend(frameon=False)
    save(fig, folder / "functional_PCA")
    result = pd.read_csv(args.output_dir / "rv_associations.csv")
    exp = (
        result[["experimental_order", "experimental_label"]]
        .drop_duplicates()
        .sort_values("experimental_order")
        .experimental_label.tolist()
    )
    sens = (
        result[["sensory_order", "sensory_label"]]
        .drop_duplicates()
        .sort_values("sensory_order")
        .sensory_label.tolist()
    )
    heat = result.pivot(
        index="experimental_label",
        columns="sensory_label",
        values="direction_encoded_rv",
    ).loc[exp, sens]
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(heat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(sens)), sens, rotation=90)
    ax.set_yticks(range(len(exp)), exp)
    q = result.pivot(
        index="experimental_label", columns="sensory_label", values="rv_q_BH_156"
    ).loc[exp, sens]
    for i in range(len(exp)):
        for j in range(len(sens)):
            ax.text(
                j,
                i,
                f"{heat.iloc[i,j]:.2f}" + ("*" if q.iloc[i, j] < 0.05 else ""),
                ha="center",
                va="center",
                fontsize=6,
            )
    fig.colorbar(im, ax=ax, label="sign(r) × RV", fraction=0.04, pad=0.02)
    fig.subplots_adjust(left=0.24, right=0.89, bottom=0.27, top=0.97)
    save(fig, folder / "direction_encoded_RV")
    tree = pd.read_csv(args.output_dir / "cluster_linkage.csv")
    fingerprint = pd.read_csv(args.output_dir / "rv_fingerprints.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    dendrogram(
        tree[["left_node", "right_node", "distance", "n_members"]].to_numpy(),
        labels=fingerprint.index.tolist(),
        ax=ax,
        leaf_rotation=90,
    )
    ax.set_ylabel("Correlation distance")
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.42, top=0.96)
    save(fig, folder / "sensory_clustering")


if __name__ == "__main__":
    main()
