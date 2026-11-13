"""Summarize the supplied LC-MS-derived structural indices and phenolic assays."""

import pandas as pd
from common import arguments


def summarize(data):
    columns = [c for c in data.select_dtypes("number") if c != "replicate"]
    out = (
        data.groupby("Fraction", observed=True)[columns]
        .agg(["mean", "std"])
        .reset_index()
    )
    out.columns = ["_".join(c).strip("_") for c in out.columns]
    return out


def main():
    args = arguments(__doc__)
    for input_name, output_name in [
        ("structure_replicates.csv", "structure_summary.csv"),
        ("total_phenolics_replicates.csv", "total_phenolics_summary.csv"),
    ]:
        summarize(pd.read_csv(args.data_dir / input_name)).to_csv(
            args.output_dir / output_name, index=False
        )


if __name__ == "__main__":
    main()
