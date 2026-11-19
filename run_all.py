"""Run analyses from public inputs and, optionally, confidential source records."""

from pathlib import Path
import argparse
import os
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=root / "results")
    parser.add_argument("--restricted-dir", type=Path)
    parser.add_argument("--restricted-output-dir", type=Path)
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def run(name, extra=(), out=None):
        subprocess.run(
            [
                sys.executable,
                str(root / "code" / name),
                "--output-dir",
                str(out or args.output_dir),
                *map(str, extra),
            ],
            check=True,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )

    for name in [
        "03_functional_pca.py",
        "04_rv_exact_permutation.py",
        "05_hierarchical_clustering.py",
        "07_qcmd_recovery_statistics.py",
        "08_physicochemical_statistics.py",
        "12_dls_distributions.py",
        "13_structure_summary.py",
        "15_bulk_endpoint_summaries.py",
    ]:
        run(name)
    if args.plots:
        run("09_plot_public_results.py")
        run("16_plot_endpoint_summaries.py")
    if args.restricted_dir:
        private_out = (
            args.restricted_output_dir
            or args.restricted_dir.resolve().parent / "results"
        )
        extra = ["--restricted-dir", args.restricted_dir]
        for name in [
            "01_ti_processing.py",
            "02_sensory_statistics.py",
            "06_qcmd_processing.py",
            "11_fpca_individual_curves.py",
            "14_score_metric_correlations.py",
        ]:
            run(name, extra, private_out)
        if args.plots:
            run("10_plot_restricted_qcmd.py", extra, private_out)
    print("Analyses completed.")


if __name__ == "__main__":
    main()
