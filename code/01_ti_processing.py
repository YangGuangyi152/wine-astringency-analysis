"""Extract the reported TI endpoints and subquality ratings from 78 assessor-fraction records."""

import json
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from common import arguments, ORDER

SUBQUALITIES = ["dry", "rough", "shrink", "gripping", "coating"]


def main():
    args = arguments(__doc__, restricted=True)
    records = [
        json.loads(line)
        for line in (args.restricted_dir / "sensory/records.jsonl")
        .read_text()
        .splitlines()
    ]
    identities = [(r["assessor_id"], r["component"]) for r in records]
    if len(records) != 78 or len(set(identities)) != 78:
        raise ValueError(
            "Expected one TI record for each of 13 assessors and six fractions."
        )
    rows = []
    curves = []
    for record in sorted(
        records, key=lambda r: (r["assessor_id"], ORDER.index(r["component"]))
    ):
        time = np.asarray(record["time_s"], dtype=float)
        raw = np.asarray(record["intensity"], dtype=float)
        dt = float(np.median(np.diff(time)))
        smooth = gaussian_filter1d(raw, sigma=1.0 / dt, mode="nearest")
        peak = int(np.argmax(smooth))
        im = float(smooth[peak])
        tm = float(time[peak])
        total = float(np.trapezoid(smooth, time))
        mask = (time >= 60) & (time <= 120)
        late = float(np.trapezoid(smooth[mask], time[mask]))
        i60 = float(np.interp(60, time, smooth))
        i90 = float(np.interp(90, time, smooth))
        half = np.where(smooth[peak:] <= 0.5 * im)[0]
        row = dict(
            user_id=record["assessor_id"],
            component=record["component"],
            Im_smooth=im,
            Tm_smooth_s=tm,
            TI_AUC_smooth=total,
            AUC_60_120=late,
            late_AUC_fraction_60_120=late / total if total > 0 else np.nan,
            I60=i60,
            I90=i90,
            half_life_after_peak_s=(
                float(time[peak + half[0]] - tm) if len(half) else np.nan
            ),
            duration_above_5_s=float(np.sum(smooth > 5) * dt),
            residual_I60_over_Im=i60 / im if im > 0 else np.nan,
        )
        row.update(
            {key + "_level": float(record["ratings"][key]) for key in SUBQUALITIES}
        )
        rows.append(row)
        curves.extend(
            dict(
                user_id=record["assessor_id"],
                component=record["component"],
                time_s=float(t),
                intensity_raw=float(y),
                intensity_smooth=float(sm),
            )
            for t, y, sm in zip(time, raw, smooth)
        )
    panel = pd.DataFrame(rows)
    panel.to_csv(args.output_dir / "participant_metrics.csv", index=False)
    long = pd.DataFrame(curves)
    long.to_csv(args.output_dir / "individual_ti_curves.csv", index=False)
    metrics = [col for col in panel if col not in ["user_id", "component"]]
    summary = panel.groupby("component")[metrics].agg(["mean", "std"]).reindex(ORDER)
    summary.columns = [
        name + "_" + ("sd" if kind == "std" else kind) for name, kind in summary.columns
    ]
    summary = summary.reset_index()
    summary.insert(1, "n_panelists", 13)
    summary.to_csv(args.output_dir / "sensory_component_summary.csv", index=False)
    long.groupby(["component", "time_s"]).intensity_smooth.agg(
        ["mean", "std", "count"]
    ).reset_index().to_csv(args.output_dir / "sensory_mean_curves.csv", index=False)


if __name__ == "__main__":
    main()
