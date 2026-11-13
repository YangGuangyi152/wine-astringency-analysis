from pathlib import Path
import argparse
import importlib.util
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ORDER = [f"F{i}" for i in range(1, 7)]


def arguments(description, restricted=False):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/public")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--restricted-dir", type=Path, required=restricted)
    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = (
            args.restricted_dir.resolve().parent / "results"
            if restricted
            else ROOT / "results"
        )
    if restricted and args.output_dir.resolve().is_relative_to(ROOT):
        parser.error("Restricted outputs must be stored outside the public repository.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    return args


def load_module(filename):
    spec = importlib.util.spec_from_file_location(
        Path(filename).stem, ROOT / "code" / filename
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    result = np.empty_like(p)
    result[order] = np.clip(adjusted, 0, 1)
    return result
