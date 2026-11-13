"""Render representative and supplementary QCM-D traces from restricted records."""

from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common import arguments
from plot_support import style, save, COLORS

COMPONENTS = ["F2", "F5", "F6"]


def smooth_series(frame, column, window_min=1 / 60):
    dt = frame.time_min.diff().median()
    points = max(5, int(round(window_min / dt)))
    if points % 2 == 0:
        points += 1
    return (
        frame[column]
        .rolling(points, center=True, min_periods=max(2, points // 4))
        .mean()
    )


def clean_axes(ax, stages):
    ax.set_xlim(stages.start_min.min(), stages.end_min.max())


def add_stage_background(ax, stages):
    for i, row in enumerate(stages.itertuples()):
        if i % 2 == 0:
            ax.axvspan(
                row.start_min, row.end_min, color="#F4F4F4", alpha=0.78, zorder=0
            )
        label = {
            "S1 rinse": "S1",
            "Tannin exposure / S2": "S2",
            "S3a rinse": "S3a",
            "S3b final rinse": "S3b",
        }[row.stage]
        ax.text(
            (row.start_min + row.end_min) / 2,
            0.985,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=10.5,
            fontweight="bold",
            color="#555555",
        )


def transition_onset(t, derivative, start, end, sign):
    indices = np.flatnonzero((t >= start) & (t <= end))
    assert len(indices) > 5
    signal = derivative * sign
    peak = indices[np.argmax(signal[indices])]
    threshold = signal[peak] * 0.12
    candidates = np.flatnonzero(
        (t >= max(start, t[peak] - 2)) & (t <= t[peak]) & (signal >= threshold)
    )
    assert len(candidates)
    return float(t[candidates[0]])


def estimate_display_boundaries(data, run, source):
    rows = []
    for component in COMPONENTS:
        frame = data[data.component == component].copy()
        t = frame.time_min.to_numpy()
        smooth = source.smooth_series(frame, "delta_f_Hz", window_min=0.5).to_numpy()
        derivative = np.gradient(smooth, t)
        s2 = transition_onset(t, derivative, 45, 70, -1)
        s3a = transition_onset(t, derivative, s2 + 12, s2 + 28, 1)
        s3b = transition_onset(t, derivative, s3a + 8, min(t[-1] - 4, s3a + 35), 1)
        assert 40 < s2 < s3a < s3b < t[-1]
        rows.append(
            {
                "run": run,
                "component": component,
                "S2_estimated_min": s2,
                "S3a_estimated_min": s3a,
                "S3b_estimated_min": s3b,
                "recording_end_min": t[-1],
            }
        )
    result = pd.DataFrame(rows)
    result["time_shift_min"] = 0.0
    medians = result[
        ["S2_estimated_min", "S3a_estimated_min", "S3b_estimated_min"]
    ].median()
    edges = [45, *medians.tolist(), 145]
    stages = pd.DataFrame(
        {
            "stage": [
                "S1 rinse",
                "Tannin exposure / S2",
                "S3a rinse",
                "S3b final rinse",
            ],
            "start_min": edges[:-1],
            "end_min": edges[1:],
        }
    )
    return (result, stages)


def draw(ax, data, metric, stages, source):
    column = "delta_f_Hz" if metric == "delta_f" else "delta_D_ppm"
    for component in COMPONENTS:
        frame = data[data.component == component]
        values = source.smooth_series(frame, column)
        ax.plot(
            frame.time_min,
            values,
            color=source.COLORS[component],
            lw=2.15,
            label=component,
        )
    ax.set_ylabel(source.METRICS[metric]["ylabel"])
    ax.set_ylim((-225, 5) if metric == "delta_f" else (-3, 50))
    source.clean_axes(ax, stages)
    ax.set_xlabel("Time (min)", fontweight="bold")
    source.add_stage_background(ax, stages)
    ax.set_xticks(np.arange(50, 141, 10))


def main():
    global COMPONENTS
    args = arguments(__doc__, restricted=True)
    style()
    plt.rcParams.update(
        {
            "axes.labelsize": 13,
            "axes.labelweight": "bold",
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
        }
    )
    source = SimpleNamespace(
        smooth_series=smooth_series,
        clean_axes=clean_axes,
        add_stage_background=add_stage_background,
        COLORS=COLORS,
        METRICS={
            "delta_f": {"ylabel": r"$\Delta f$ (Hz)"},
            "delta_d": {"ylabel": r"$\Delta D$ (ppm)"},
        },
    )
    for system in ["BSM", "saliva"]:
        COMPONENTS = (
            [f"F{i}" for i in range(1, 7)] if system == "BSM" else ["F2", "F5", "F6"]
        )
        recordings = {}
        stage_tables = {}
        bounds = []
        for run in ["1", "2", "3"]:
            frames = []
            for c in COMPONENTS:
                d = pd.read_csv(
                    args.restricted_dir / "qcmd" / system / f"{run}_{c}.csv"
                )
                d["time_min"] = d.time_s / 60
                d["component"] = c
                frames.append(d)
            data = pd.concat(frames, ignore_index=True)
            boundary, stages = estimate_display_boundaries(data, run, source)
            if system == "BSM" and run != "1":
                shift = boundary.set_index("component").S2_estimated_min
                data["time_min"] = data.time_min - data.component.map(shift)
                edges = [
                    -12,
                    0,
                    float(
                        (
                            boundary.S3a_estimated_min - boundary.S2_estimated_min
                        ).median()
                    ),
                    float(
                        (
                            boundary.S3b_estimated_min - boundary.S2_estimated_min
                        ).median()
                    ),
                    75,
                ]
                stages = pd.DataFrame(
                    {
                        "stage": [
                            "S1 rinse",
                            "Tannin exposure / S2",
                            "S3a rinse",
                            "S3b final rinse",
                        ],
                        "start_min": edges[:-1],
                        "end_min": edges[1:],
                    }
                )
                boundary["time_shift_min"] = -boundary.S2_estimated_min
            bounds.append(boundary)
            recordings[run] = data
            stage_tables[run] = stages
        for runs, stem in [
            (["1"], system + "_representative_experiment"),
            (["2", "3"], system + "_additional_experiments"),
        ]:
            fig, axes = plt.subplots(
                len(runs), 2, figsize=(14.2, 4.15 * len(runs)), squeeze=False
            )
            fig.subplots_adjust(
                left=0.08,
                right=0.98,
                bottom=0.08 if len(runs) == 2 else 0.18,
                top=0.85,
                wspace=0.23,
                hspace=0.31,
            )
            for row, run in enumerate(runs):
                for col, metric in enumerate(["delta_f", "delta_d"]):
                    ax = axes[row, col]
                    draw(ax, recordings[run], metric, stage_tables[run], source)
                    if metric == "delta_d" and system == "BSM":
                        ax.set_ylim(-3, 130)
                    if system == "BSM" and run != "1":
                        ax.set_xlabel(
                            "Time relative to S2 onset (min)", fontweight="bold"
                        )
                        ax.set_xticks(np.arange(-10, 71, 10))
                    ax.set_xlim(
                        stage_tables[run].start_min.min(),
                        stage_tables[run].end_min.max(),
                    )
                    ax.text(
                        -0.13,
                        1.04,
                        "ABCD"[row * 2 + col],
                        transform=ax.transAxes,
                        fontweight="bold",
                        fontsize=16,
                        ha="left",
                        va="bottom",
                    )
            fig.legend(
                *axes[0, 0].get_legend_handles_labels(),
                loc="upper center",
                bbox_to_anchor=(0.53, 0.995),
                ncol=len(COMPONENTS),
                frameon=False,
            )
            save(fig, args.output_dir / stem)
        pd.concat(bounds, ignore_index=True).to_csv(
            args.output_dir / f"{system}_trace_display_boundaries.csv", index=False
        )


if __name__ == "__main__":
    main()
