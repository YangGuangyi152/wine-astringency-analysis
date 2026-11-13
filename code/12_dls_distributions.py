"""Prepare DLS distribution curves with the original log-grid display smoothing."""

import numpy as np
import pandas as pd
from common import arguments

MIN_PEAK_DIAMETER = 90
MAX_PEAK_DIAMETER = 10200


def gaussian_kernel(size=21, sigma=2.2):
    grid = np.arange(size) - (size - 1) / 2
    kernel = np.exp(-0.5 * (grid / sigma) ** 2)
    return kernel / kernel.sum()


def gaussian_smoothed_curve(x, y, x_min, x_max, n_points=650):
    if len(x) < 4:
        return (x, y)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    y = np.where((x >= MIN_PEAK_DIAMETER) & (x <= MAX_PEAK_DIAMETER), y, 0)
    positive = y > 0
    if positive.any():
        first = np.argmax(positive)
        last = len(positive) - np.argmax(positive[::-1])
        support_x = x[first:last]
        support_y = y[first:last]
    else:
        support_x = x
        support_y = y
    log_x = np.log10(x)
    dense_log_x = np.linspace(np.log10(x_min), np.log10(x_max), n_points)
    dense_y = np.interp(dense_log_x, log_x, y, left=0, right=0)
    kernel = gaussian_kernel()
    padded = np.pad(
        dense_y,
        (len(kernel) // 2, len(kernel) // 2),
        mode="constant",
        constant_values=0,
    )
    smooth_y = np.convolve(padded, kernel, mode="valid")
    if smooth_y.max() > 0 and support_y.max() > 0:
        smooth_y = smooth_y / smooth_y.max() * support_y.max()
    smooth_y[0] = 0
    smooth_y[-1] = 0
    return (10**dense_log_x, smooth_y)


def group_axis_limits(curves, group):
    xs = []
    for prefix in [group, f"{group}-3h"]:
        for repeat_no in [1, 2, 3]:
            x, y = curves[prefix, repeat_no]
            positive_x = x[
                (y > 0) & (x >= MIN_PEAK_DIAMETER) & (x <= MAX_PEAK_DIAMETER)
            ]
            if len(positive_x):
                xs.extend([positive_x.min(), positive_x.max()])
    if not xs:
        return (MIN_PEAK_DIAMETER, MAX_PEAK_DIAMETER)
    min_power = int(np.floor(np.log10(min(xs))))
    max_power = int(np.ceil(np.log10(max(xs))))
    xmin = max(MIN_PEAK_DIAMETER, 10**min_power)
    xmax = min(MAX_PEAK_DIAMETER, 10**max_power)
    return (xmin, xmax)


def main():
    args = arguments(__doc__)
    table = pd.read_csv(args.data_dir / "dls_size_distributions.csv")
    curves = {
        (group, int(rep)): (g.diameter_nm.to_numpy(), g.intensity.to_numpy())
        for (group, rep), g in table.groupby(["sample_group", "replicate"])
    }
    rows = []
    for group in [f"BSM-F{i}" for i in range(1, 7)]:
        lo, hi = group_axis_limits(curves, group)
        for sample in [group, group + "-3h"]:
            for rep in [1, 2, 3]:
                x, y = curves[(sample, rep)]
                xx, yy = gaussian_smoothed_curve(x, y, lo, hi)
                rows.extend(
                    dict(sample_group=sample, replicate=rep, diameter_nm=a, intensity=b)
                    for a, b in zip(xx, yy)
                )
    pd.DataFrame(rows).to_csv(args.output_dir / "dls_display_curves.csv", index=False)


if __name__ == "__main__":
    main()
