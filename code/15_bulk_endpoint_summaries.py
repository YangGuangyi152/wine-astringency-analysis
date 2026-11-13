"""Summarize DLS and electrophoretic endpoints without changing missing measurements."""

import re
import pandas as pd
from common import arguments


def main():
    args = arguments(__doc__)
    d = pd.read_csv(args.data_dir / "dls_summary_replicates.csv")
    d["sample_group"] = (
        d["Sample ID"].str.replace(r"\s*-\s*\d+\s*$", "", regex=True).str.strip()
    )
    columns = [c for c in d.select_dtypes("number")]
    out = (
        d.groupby("sample_group", sort=True)[columns]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    out.columns = ["_".join(c).strip("_") for c in out.columns]
    out.to_csv(args.output_dir / "dls_endpoint_summary.csv", index=False)
    z = pd.read_csv(args.data_dir / "zeta_potential_replicates.csv")
    columns = [
        "Zeta Potential (mV)",
        "Mobility (μ/s)/(V/cm)",
        "Conductance (μS)",
        "RMS Residual",
    ]
    out = z.groupby("group")[columns].agg(["mean", "std", "count"]).reset_index()
    out.columns = ["_".join(col).strip("_") for col in out.columns]
    out.to_csv(args.output_dir / "electrokinetic_summary.csv", index=False)


if __name__ == "__main__":
    main()
