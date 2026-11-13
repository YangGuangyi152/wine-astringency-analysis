"""Summarize processed per-experiment recoveries and extract saliva stage means."""

import pandas as pd
from common import arguments


def main():
    args = arguments(__doc__, restricted=True)
    summaries = []
    for system in ["BSM", "saliva"]:
        folder = args.restricted_dir / "qcmd" / system
        recovery = pd.read_csv(folder / "recovery.csv").assign(film_system=system)
        recovery.to_csv(args.output_dir / f"{system}_recovery.csv", index=False)
        summaries.append(
            recovery.groupby(["film_system", "component", "metric", "stage"])
            .recovery_percent.agg(["mean", "std", "count"])
            .reset_index()
        )
        if system == "saliva":
            windows = pd.read_csv(folder / "stage_windows.csv")
            rows = []
            for w in windows.itertuples():
                frame = pd.read_csv(folder / f"{int(w.run)}_{w.component}.csv")
                time = frame.time_s / 60
                selected = frame[
                    (time >= w.window_start_min) & (time < w.stage_end_min)
                ]
                for metric, col in [
                    ("delta_f", "delta_f_Hz"),
                    ("delta_d", "delta_D_ppm"),
                ]:
                    rows.append(
                        dict(
                            run=int(w.run),
                            component=w.component,
                            metric=metric,
                            stage=w.stage,
                            mean=selected[col].mean(),
                        )
                    )
            pd.DataFrame(rows).to_csv(
                args.output_dir / "saliva_stage_means_from_traces.csv", index=False
            )
    pd.concat(summaries, ignore_index=True).to_csv(
        args.output_dir / "qcmd_recovery_summary_from_private.csv", index=False
    )


if __name__ == "__main__":
    main()
