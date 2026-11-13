"""Reproduce the centered full-curve PCA from aggregate second moments."""

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.stats import f
from common import arguments, ORDER


def main():
    args = arguments(__doc__)
    data = np.load(args.data_dir / "fpca_aggregate_moments.npz")
    cov = data["covariance"]
    vals, vec = eigh(cov, subset_by_index=[len(cov) - 2, len(cov) - 1])
    vals = vals[::-1]
    loadings = vec[:, ::-1].T
    for row in loadings:
        if row[np.argmax(np.abs(row))] < 0:
            row *= -1
    variance = pd.DataFrame(
        {
            "PC": ["PC1", "PC2"],
            "explained_variance": vals,
            "explained_variance_ratio": vals / np.trace(cov),
        }
    )
    variance.to_csv(args.output_dir / "fpca_variance.csv", index=False)
    pd.DataFrame(
        {
            "time_s": data["time_s"],
            "PC1_loading": loadings[0],
            "PC2_loading": loadings[1],
        }
    ).to_csv(args.output_dir / "fpca_loadings.csv", index=False)
    scores = data["component_centered_mean"] @ loadings.T
    moments = pd.read_csv(args.data_dir / "fpca_group_score_covariance.csv").set_index(
        "component"
    )
    rows = []
    for i, c in enumerate(ORDER):
        m = moments.loc[c]
        n = int(m["n"])
        covariance = np.array(
            [[m.cov_PC1_PC1, m.cov_PC1_PC2], [m.cov_PC1_PC2, m.cov_PC2_PC2]]
        )
        factor = 2 * (n - 1) / (n - 2) * f.ppf(0.95, 2, n - 2)
        eigen, vectors = np.linalg.eigh(covariance / n)
        order = np.argsort(eigen)[::-1]
        eigen = eigen[order]
        vectors = vectors[:, order]
        semi = np.sqrt(factor * np.maximum(eigen, 0))
        rows.append(
            dict(
                component=c,
                n=n,
                PC1_mean=scores[i, 0],
                PC2_mean=scores[i, 1],
                cov_PC1_PC1=covariance[0, 0],
                cov_PC1_PC2=covariance[0, 1],
                cov_PC2_PC2=covariance[1, 1],
                hotelling_constant_95=factor,
                ellipse_width_95=2 * semi[0],
                ellipse_height_95=2 * semi[1],
                ellipse_angle_deg=np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0])),
                PC1_projected_radius_95=np.sqrt(factor * covariance[0, 0] / n),
                PC2_projected_radius_95=np.sqrt(factor * covariance[1, 1] / n),
            )
        )
    pd.DataFrame(rows).to_csv(
        args.output_dir / "fpca_centroids_and_ellipses.csv", index=False
    )


if __name__ == "__main__":
    main()
