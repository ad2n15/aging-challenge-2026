#!/usr/bin/env python3
"""
Verify that public competition .h5ad files do not leak ground-truth age for the test split.

Checks these objects relative to the chosen data root (see --data-dir). For each
logical file, uses ``*_public.h5ad`` if present (output of ``strip_test_age_h5ad.py``),
otherwise the canonical name (e.g. ``combined.h5ad``).

Default data root: prefers ``data/``, then ``data_prep/output_public/`` (stripped
``*_public.h5ad``), then ``data_prep/output/``, else ``data/``.

For every row with obs['_split'] == 'test', column 'age' must be entirely missing/NaN
(if the column exists). Any finite age value on test rows is treated as a leak.

Requires: scanpy, pandas, numpy (same as teaching notebooks).

Usage:
  python scripts/check_test_age_withheld.py
  python scripts/check_test_age_withheld.py --data-dir /path/to/data

Exit codes:
  0 — all present files pass (no test-age leak)
  1 — missing file, missing _split, or leak detected
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def _default_data_dir() -> Path:
    root = Path(__file__).resolve().parents[1]
    data = root / "data"
    prep_out = root / "data_prep" / "output"
    out_pub = root / "data_prep" / "output_public"
    if (data / "combined.h5ad").is_file() or (data / "combined_public.h5ad").is_file():
        return data
    if (out_pub / "combined_public.h5ad").is_file() or (out_pub / "combined.h5ad").is_file():
        return out_pub
    if (prep_out / "combined.h5ad").is_file() or (prep_out / "combined_public.h5ad").is_file():
        return prep_out
    return data


def _h5ad_path(data_dir: Path, rel: Path, name_suffix: str) -> Path:
    """Prefer ``{stem}{suffix}.h5ad`` if that file exists, else ``rel``."""
    tagged = rel.with_name(f"{rel.stem}{name_suffix}{rel.suffix}")
    p_tag = data_dir / tagged
    p_std = data_dir / rel
    if p_tag.is_file():
        return p_tag
    return p_std


def _test_mask(obs: pd.DataFrame) -> pd.Series:
    if "_split" not in obs.columns:
        raise ValueError("obs has no '_split' column")
    return obs["_split"].astype(str).str.strip().str.lower() == "test"


def _ages_for_test(obs: pd.DataFrame, test_mask: pd.Series) -> pd.Series | None:
    for col in ("age", "Age", "AGE"):
        if col in obs.columns:
            return pd.to_numeric(obs.loc[test_mask, col], errors="coerce")
    return None


def check_h5ad(path: Path, label: str) -> tuple[bool, list[str]]:
    """Return (ok, lines of report)."""
    lines: list[str] = []
    if not path.is_file():
        lines.append(f"{label}: ERROR — file not found: {path}")
        return False, lines

    adata = sc.read_h5ad(path, backed="r")
    try:
        obs = adata.obs
        n_rows = obs.shape[0]
        try:
            tm = _test_mask(obs)
        except ValueError as e:
            lines.append(f"{label}: ERROR — {e}")
            return False, lines
        n_test = int(tm.sum())
        lines.append(f"{label}: {path.name} — {n_rows} rows, {n_test} with _split=='test'")

        if n_test == 0:
            lines.append(f"  WARNING: no test split rows; check file version.")
            return True, lines

        ages = _ages_for_test(obs, tm)
        if ages is None:
            lines.append("  OK — no 'age' column (nothing to leak on test).")
            return True, lines

        vals = ages.to_numpy(dtype=float)
        finite = np.isfinite(vals)
        n_leak = int(finite.sum())
        if n_leak > 0:
            lines.append(f"  FAIL — {n_leak} test rows have finite age values (leak).")
            return False, lines

        lines.append("  OK — 'age' column exists but all test values are NaN.")
        return True, lines
    finally:
        if getattr(adata, "isbacked", False) and adata.file is not None:
            adata.file.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Folder with combined*.h5ad and pseudobulk/ (default: auto data/, output_public, or data_prep/output/)",
    )
    parser.add_argument(
        "--name-suffix",
        type=str,
        default="_public",
        help="Prefer files named {stem}{suffix}.h5ad when present (default: _public). Use '' to only check canonical names.",
    )
    args = parser.parse_args()
    data_dir = args.data_dir if args.data_dir is not None else _default_data_dir()
    sfx = args.name_suffix

    rels = [
        Path("combined.h5ad"),
        Path("pseudobulk") / "combined_pseudobulk_combined.h5ad",
        Path("pseudobulk") / "combined_pseudobulk_donor_aggregated.h5ad",
    ]
    files: list[tuple[str, Path]] = [
        ("Cell-level combined", _h5ad_path(data_dir, rels[0], sfx)),
        ("Pseudobulk combined (donor × cell type)", _h5ad_path(data_dir, rels[1], sfx)),
        ("Pseudobulk donor-aggregated", _h5ad_path(data_dir, rels[2], sfx)),
    ]

    print(f"Data directory: {data_dir.resolve()}\n")
    all_ok = True
    for label, p in files:
        ok, lines = check_h5ad(p, label)
        for line in lines:
            print(line)
        print()
        if not ok:
            all_ok = False

    if all_ok:
        print("Summary: PASS — no test-split age leak detected in checked files.")
        return 0
    print("Summary: FAIL — fix data or paths before publishing.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
