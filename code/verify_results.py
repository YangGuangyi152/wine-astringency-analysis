"""Verify manuscript numerical outputs against frozen reference tables."""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from common import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/public")
    parser.add_argument("--restricted-dir", type=Path)
    parser.add_argument("--restricted-output-dir", type=Path)
    args = parser.parse_args()
    checks = []
    ref = ROOT / "reference_results"
    out = args.output_dir

    def check(label, expected, observed, keys, columns=None):
        a = pd.read_csv(expected).sort_values(keys).reset_index(drop=True)
        b = pd.read_csv(observed).sort_values(keys).reset_index(drop=True)
        assert len(a) == len(b), label
        for key in keys:
            assert a[key].astype(str).tolist() == b[key].astype(str).tolist(), (
                label + "/" + key
            )
        numeric = columns or [
            col for col in a.select_dtypes("number") if col in b and col not in keys
        ]
        errors = []
        for col in numeric:
            np.testing.assert_allclose(
                a[col],
                b[col],
                rtol=1e-10,
                atol=1e-12,
                equal_nan=True,
                err_msg=label + "/" + col,
            )
            delta = np.abs(a[col].to_numpy(float) - b[col].to_numpy(float))
            errors.extend(delta[np.isfinite(delta)])
        checks.append(
            dict(
                check=label,
                status="PASS",
                rows=len(a),
                numeric_columns=len(numeric),
                max_absolute_error=max(errors, default=0.0),
            )
        )

    for name, keys in [
        ("fpca_variance.csv", ["PC"]),
        ("fpca_loadings.csv", ["time_s"]),
        ("fpca_centroids_and_ellipses.csv", ["component"]),
    ]:
        check("Public " + name, ref / name, out / name, keys)
    check(
        "156 exact-permutation RV associations",
        ref / "rv_associations.csv",
        out / "rv_associations.csv",
        ["experimental_column", "sensory_column"],
    )
    check(
        "Clustering linkage",
        ref / "cluster_linkage.csv",
        out / "cluster_linkage.csv",
        ["merge_step"],
    )
    check(
        "Average-linkage silhouette screen",
        ref / "clustering_parameter_screen.csv",
        out / "clustering_parameter_screen.csv",
        ["data_basis", "linkage", "k"],
        ["silhouette"],
    )
    a = pd.read_csv(ref / "cluster_assignments.csv").sort_values("sensory_indicator")
    b = pd.read_csv(out / "cluster_assignments.csv").sort_values("sensory_indicator")
    assert a.empirical_module.tolist() == b.empirical_module.tolist()
    checks.append(
        dict(check="Four empirical module memberships", status="PASS", rows=len(a))
    )
    for system in ["BSM", "saliva"]:
        for stem, keys in [
            ("recovery_summary", ["metric", "stage", "component"]),
            ("recovery_anova", ["metric", "stage"]),
            ("recovery_tukey", ["metric", "stage", "group_a", "group_b"]),
        ]:
            check(
                system + " " + stem,
                ref / f"{system}_{stem}.csv",
                out / f"{system}_{stem}.csv",
                keys,
            )
        a = pd.read_csv(ref / f"{system}_recovery_summary.csv").sort_values(
            ["metric", "stage", "component"]
        )
        b = pd.read_csv(out / f"{system}_recovery_summary.csv").sort_values(
            ["metric", "stage", "component"]
        )
        assert a.letters.tolist() == b.letters.tolist()
    for name in ["structure_summary.csv", "total_phenolics_summary.csv"]:
        check(name, args.data_dir / name, out / name, ["Fraction"])
    check(
        "Bulk ANOVA/Tukey comparisons",
        ref / "bulk_ordinary_ANOVA_Tukey_pairs.csv",
        out / "bulk_pairwise.csv",
        ["endpoint", "group_a", "group_b"],
    )
    expected = (
        pd.read_csv(ref / "bulk_assumption_diagnostics_and_omnibus.csv")
        .set_index("endpoint")
        .sort_index()
    )
    diagnostic = (
        pd.read_csv(out / "bulk_diagnostics.csv").set_index("endpoint").sort_index()
    )
    omnibus = pd.read_csv(out / "bulk_omnibus.csv").set_index("endpoint").sort_index()
    for old, new in [
        ("shapiro_W", "shapiro_W"),
        ("shapiro_p", "shapiro_p"),
        ("levene_median_F", "levene_F"),
        ("levene_median_p", "levene_p"),
    ]:
        np.testing.assert_allclose(
            expected[old], diagnostic[new], rtol=1e-10, atol=1e-12
        )
    for old, new in [("ordinary_F", "F"), ("ordinary_p", "p_value")]:
        np.testing.assert_allclose(expected[old], omnibus[new], rtol=1e-10, atol=1e-12)
    checks.append(
        dict(
            check="Bulk assumption diagnostics and omnibus tests",
            status="PASS",
            rows=len(expected),
        )
    )
    variance = pd.read_csv(out / "fpca_variance.csv")
    assert [f"{x * 100:.1f}" for x in variance.explained_variance_ratio] == [
        "81.1",
        "14.9",
    ]
    silhouette = (
        pd.read_csv(out / "clustering_parameter_screen.csv")
        .query("k==4")
        .silhouette.iloc[0]
    )
    assert f"{silhouette:.3f}" == "0.790"
    checks.append(
        dict(check="Manuscript-rounded fPCA and silhouette values", status="PASS")
    )
    if args.restricted_dir:
        private = (
            args.restricted_output_dir
            or args.restricted_dir.resolve().parent / "results"
        )
        if (args.restricted_dir / "sensory/participant_metrics_reference.csv").exists():
            check(
                "Reported individual sensory endpoints",
                args.restricted_dir / "sensory/participant_metrics_reference.csv",
                private / "participant_metrics.csv",
                ["user_id", "component"],
            )
        check(
            "Reported paired sensory comparisons",
            ref / "sensory_pairwise_reference.csv",
            private / "sensory_pairwise.csv",
            ["metric", "comparison"],
        )
        check(
            "Friedman and Kendall W",
            ref / "sensory_friedman_reference.csv",
            private / "sensory_friedman.csv",
            ["metric"],
        )
        for name, keys in [
            ("fpca_variance.csv", ["PC"]),
            ("fpca_loadings.csv", ["time_s"]),
            ("fpca_centroids_and_ellipses.csv", ["component"]),
        ]:
            check("Individual-data " + name, ref / name, private / name, keys)
        if (args.restricted_dir / "sensory/fpca_scores_reference.csv").exists():
            check(
                "Individual fPC scores",
                args.restricted_dir / "sensory/fpca_scores_reference.csv",
                private / "fpca_individual_scores.csv",
                ["assessor_id", "component"],
            )
        check(
            "Five cited fPC1 correlations",
            ref / "score_metric_spearman_reference.csv",
            private / "score_metric_spearman_correlations.csv",
            ["PC", "metric"],
        )
        for name, keys in [
            ("sensory_component_summary.csv", ["component"]),
            ("sensory_mean_curves.csv", ["component", "time_s"]),
        ]:
            check("Recomputed " + name, args.data_dir / name, private / name, keys)
        expected = np.load(args.data_dir / "fpca_aggregate_moments.npz")
        computed = np.load(private / "fpca_aggregate_moments_recomputed.npz")
        for key in expected.files:
            np.testing.assert_allclose(
                expected[key], computed[key], rtol=1e-10, atol=1e-12
            )
        checks.append(
            dict(
                check="PCA aggregate moments from individual curves",
                status="PASS",
                rows=78,
            )
        )
        for system in ["BSM", "saliva"]:
            if (args.restricted_dir / f"qcmd/{system}/recovery.csv").exists():
                check(
                    system + " per-experiment recovery input",
                    args.restricted_dir / f"qcmd/{system}/recovery.csv",
                    private / f"{system}_recovery.csv",
                    ["run", "component", "metric", "stage"],
                    ["recovery_percent"],
                )
        if (private / "qcmd_recovery_summary_from_private.csv").exists():
            check(
                "Recovery summaries from three independent experiments",
                args.data_dir / "qcmd_recovery_summary.csv",
                private / "qcmd_recovery_summary_from_private.csv",
                ["film_system", "component", "metric", "stage"],
            )
        if (args.restricted_dir / "qcmd/saliva/stage_means.csv").exists() and (
            private / "saliva_stage_means_from_traces.csv"
        ).exists():
            check(
                "Saliva stage means from traces",
                args.restricted_dir / "qcmd/saliva/stage_means.csv",
                private / "saliva_stage_means_from_traces.csv",
                ["run", "component", "metric", "stage"],
                ["mean"],
            )
    pd.DataFrame(checks).to_csv(out / "verification_report.csv", index=False)
    (out / "verification_report.json").write_text(json.dumps(checks, indent=2))
    print(f"PASS: {len(checks)} numerical verification checks.")


if __name__ == "__main__":
    main()
