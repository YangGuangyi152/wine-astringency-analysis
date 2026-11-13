"""Portable figure exports and optional local panel-alignment auditing."""

import json
import os
import importlib.util
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
    audit = os.environ.get("FIGURE_ALIGNMENT_AUDITOR")
    if audit:
        spec = importlib.util.spec_from_file_location("alignment_auditor", audit)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.require_matplotlib_panel_alignment(
            fig, json_out=str(path.with_suffix(".alignment.json")), strict=True
        )
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".svg"))
    fig.savefig(path.with_suffix(".png"), dpi=200)
    fig.savefig(
        path.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"}
    )
    plt.close(fig)
