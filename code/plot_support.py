"""Export figures in PDF, SVG, PNG and TIFF formats."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {
    "F1": "#2F3A78",
    "F2": "#6FA8C9",
    "F3": "#8F2F4B",
    "F4": "#D8A0B9",
    "F5": "#6FA56E",
    "F6": "#D39A2D",
}


def style():
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": None,
        }
    )


def save(fig, path):
    fig.canvas.draw()
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".svg"))
    fig.savefig(path.with_suffix(".png"), dpi=200)
    fig.savefig(
        path.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"}
    )
    plt.close(fig)
