"""Run the public workflow, or additionally verify reviewer-only individual data."""

from pathlib import Path
import argparse
import subprocess
import sys
import time
import json


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=root / "results")
    parser.add_argument("--restricted-dir", type=Path)
    parser.add_argument("--restricted-output-dir", type=Path)
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    timings = []

    def run(name, extra=(), out=None):
        start = time.perf_counter()
        subprocess.run(
            [
                sys.executable,
                str(root / "code" / name),
                "--output-dir",
                str(out or args.output_dir),
                *map(str, extra),
            ],
            check=True,
        )
        timings.append(dict(script=name, seconds=time.perf_counter() - start))

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
        run(
            "verify_results.py",
            [
                "--restricted-dir",
                args.restricted_dir,
                "--restricted-output-dir",
                private_out,
            ],
        )
    else:
        run("verify_results.py")
    (args.output_dir / "run_times.json").write_text(json.dumps(timings, indent=2))
    print("Completed all requested analyses and numerical verification.")


if __name__ == "__main__":
    main()
