"""Plot public structural, bulk and subquality summaries from processed inputs."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common import arguments, ORDER
from plot_support import style, save, COLORS


def main():
    args = arguments(__doc__)
    style()
    folder = args.output_dir / "figures"
    folder.mkdir(exist_ok=True)
    structure = (
        pd.read_csv(args.output_dir / "structure_summary.csv")
        .set_index("Fraction")
        .loc[ORDER]
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.3, 3.4))
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.17, top=0.96, wspace=0.4)
    for ax, (col, label) in zip(
        axes,
        [
            ("mDP", "Mean degree of polymerization"),
            ("Extension_percent", "Extension units (%)"),
            ("PD_percent", "Prodelphinidin (%)"),
        ],
    ):
        ax.bar(
            ORDER,
            structure[f"{col}_mean"],
            yerr=structure[f"{col}_std"],
            color=[COLORS[c] for c in ORDER],
            capsize=3,
        )
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0)
    save(fig, folder / "structure_summaries")
    bulk = pd.read_csv(args.output_dir / "bulk_summary.csv").query(
        "method=='ANOVA_Tukey'"
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.3, 3.7))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.2, top=0.94, wspace=0.35)
    for ax, (endpoint, label) in zip(
        axes,
        [
            ("A600", "Absorbance at 600 nm"),
            ("DLS", "Apparent hydrodynamic diameter (nm)"),
        ],
    ):
        rows = (
            bulk[bulk.endpoint == endpoint]
            .set_index("group")
            .reindex(["BSM"] + [f"BSM-{c}" for c in ORDER])
        )
        x = np.arange(len(rows))
        ax.bar(
            x,
            rows["mean"],
            yerr=rows["std"],
            color=["#B7B7B7"] + [COLORS[c] for c in ORDER],
            capsize=3,
        )
        ax.set_xticks(x, ["BSM"] + ORDER, rotation=35, ha="right")
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0)
    save(fig, folder / "bulk_endpoint_summaries")
    curves = (
        pd.read_csv(args.output_dir / "dls_display_curves.csv")
        .groupby(["sample_group", "diameter_nm"])
        .intensity.mean()
        .reset_index()
    )
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.5))
    fig.subplots_adjust(
        left=0.09, right=0.98, bottom=0.12, top=0.94, wspace=0.32, hspace=0.36
    )
    for ax, c in zip(axes.flat, ORDER):
        for label, ls in [(f"BSM-{c}", "-"), (f"BSM-{c}-3h", "--")]:
            subset = curves[curves.sample_group == label]
            if not subset.empty:
                ax.plot(
                    subset.diameter_nm,
                    subset.intensity,
                    color=COLORS[c],
                    ls=ls,
                    label="Initial" if ls == "-" else "3 h",
                )
        ax.set_xscale("log")
        ax.set(xlabel="Diameter (nm)", ylabel="Intensity", title=c)
        ax.legend(frameon=False, fontsize=8)
    save(fig, folder / "DLS_distributions")
    sensory = (
        pd.read_csv(args.data_dir / "sensory_component_summary.csv")
        .set_index("component")
        .loc[ORDER]
    )
    attributes = [
        ("dry_level", "Drying"),
        ("rough_level", "Rough"),
        ("shrink_level", "Pucker"),
        ("gripping_level", "Grippy"),
        ("coating_level", "Mouthcoating"),
    ]
    fig, axes = plt.subplots(1, 5, figsize=(12.6, 3.1))
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.16, top=0.9, wspace=0.36)
    for ax, (col, label) in zip(axes, attributes):
        ax.bar(
            ORDER,
            sensory[f"{col}_mean"],
            yerr=sensory[f"{col}_sd"],
            color=[COLORS[c] for c in ORDER],
            capsize=2,
        )
        ax.set(title=label, ylim=(0, 5.6), yticks=[1, 2, 3, 4, 5])
    axes[0].set_ylabel("Intensity category")
    save(fig, folder / "subquality_summaries")


if __name__ == "__main__":
    main()
