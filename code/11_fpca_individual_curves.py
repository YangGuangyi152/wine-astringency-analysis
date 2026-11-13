"""Full-curve PCA and joint confidence ellipses from restricted individual curves."""

import numpy as np
import pandas as pd
from scipy.linalg import svd
from scipy.stats import t, f
from common import arguments, ORDER


def confidence_tables(
    scores: pd.DataFrame, alpha: float = 0.05
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cross_rows = []
    ellipse_rows = []
    p_dim = 2
    for component in ORDER:
        values = scores.loc[scores["component"].eq(component), ["PC1", "PC2"]].to_numpy(
            float
        )
        n = len(values)
        mean = values.mean(axis=0)
        sd = values.std(axis=0, ddof=1)
        sem = sd / np.sqrt(n)
        t_critical = float(t.ppf(1 - alpha / 2, n - 1))
        margin = t_critical * sem
        cross_rows.append(
            {
                "component": component,
                "n": n,
                "PC1_mean": mean[0],
                "PC2_mean": mean[1],
                "PC1_sd": sd[0],
                "PC2_sd": sd[1],
                "PC1_sem": sem[0],
                "PC2_sem": sem[1],
                "t_critical_95": t_critical,
                "PC1_margin_95": margin[0],
                "PC2_margin_95": margin[1],
                "PC1_ci95_low": mean[0] - margin[0],
                "PC1_ci95_high": mean[0] + margin[0],
                "PC2_ci95_low": mean[1] - margin[1],
                "PC2_ci95_high": mean[1] + margin[1],
            }
        )
        covariance = np.cov(values, rowvar=False, ddof=1)
        covariance_mean = covariance / n
        hotelling_constant = float(
            p_dim * (n - 1) / (n - p_dim) * f.ppf(1 - alpha, p_dim, n - p_dim)
        )
        eigenvalues, eigenvectors = np.linalg.eigh(covariance_mean)
        eigenvalues = np.maximum(eigenvalues, 0)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]
        semi_axes = np.sqrt(hotelling_constant * eigenvalues)
        angle = float(np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])))
        ellipse_rows.append(
            {
                "component": component,
                "n": n,
                "PC1_mean": mean[0],
                "PC2_mean": mean[1],
                "cov_PC1_PC1": covariance[0, 0],
                "cov_PC1_PC2": covariance[0, 1],
                "cov_PC2_PC2": covariance[1, 1],
                "hotelling_constant_95": hotelling_constant,
                "ellipse_width_95": 2 * semi_axes[0],
                "ellipse_height_95": 2 * semi_axes[1],
                "ellipse_angle_deg": angle,
                "PC1_projected_radius_95": np.sqrt(
                    hotelling_constant * covariance_mean[0, 0]
                ),
                "PC2_projected_radius_95": np.sqrt(
                    hotelling_constant * covariance_mean[1, 1]
                ),
            }
        )
    return (pd.DataFrame(cross_rows), pd.DataFrame(ellipse_rows))


def main():
    args = arguments(__doc__, restricted=True)
    long = pd.read_csv(args.output_dir / "individual_ti_curves.csv")
    matrix = long.pivot(
        index=["component", "user_id"], columns="time_s", values="intensity_smooth"
    ).sort_index()
    x = matrix.to_numpy(float)
    centered = x - x.mean(0)
    u, s, vt = svd(centered, full_matrices=False)
    for i, row in enumerate(vt):
        if row[np.argmax(np.abs(row))] < 0:
            vt[i] *= -1
            u[:, i] *= -1
    scores = matrix.index.to_frame(index=False)
    scores["PC1"] = u[:, 0] * s[0]
    scores["PC2"] = u[:, 1] * s[1]
    scores = scores.rename(columns={"user_id": "assessor_id"})
    scores.to_csv(args.output_dir / "fpca_individual_scores.csv", index=False)
    time_labels = np.asarray([float(f"{x:.1f}") for x in matrix.columns])
    pd.DataFrame(
        {"time_s": time_labels, "PC1_loading": vt[0], "PC2_loading": vt[1]}
    ).to_csv(args.output_dir / "fpca_loadings.csv", index=False)
    pd.DataFrame(
        {
            "PC": ["PC1", "PC2"],
            "explained_variance": s[:2] ** 2 / (len(x) - 1),
            "explained_variance_ratio": s[:2] ** 2 / np.sum(s**2),
        }
    ).to_csv(args.output_dir / "fpca_variance.csv", index=False)
    _, ellipse = confidence_tables(scores)
    ellipse.to_csv(args.output_dir / "fpca_centroids_and_ellipses.csv", index=False)
    np.savez_compressed(
        args.output_dir / "fpca_aggregate_moments_recomputed.npz",
        time_s=time_labels,
        covariance=centered.T @ centered / (len(centered) - 1),
        component_centered_mean=np.vstack(
            [
                centered[matrix.index.get_level_values("component") == c].mean(0)
                for c in ORDER
            ]
        ),
        component_n=np.array([13] * 6),
        n_observations=len(centered),
    )


if __name__ == "__main__":
    main()
