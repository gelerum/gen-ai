#!/usr/bin/env python3
"""Clean the RCSB PDB per-chain features into a foldability training set.

Input is pdb_protein_features.parquet from extract_pdb_features.py: one row per
chain with a per-residue CA `bfactor` array, plus `method`, `resolution`,
`wilson_b`, `sequence` and taxonomy. We keep only well-measured X-ray chains and
normalise B-factors by the structure's Wilson B so they compare across entries.

Steps (each logged with a before/after row count):
  1. keep method == "X-RAY DIFFRACTION" (Wilson B is only defined for X-ray)
  2. keep rows with a usable wilson_b (not null and > 0)
  3. bfactor_norm = bfactor / wilson_b, then per chain:
       - bnorm_mean  = mean of bfactor_norm over residues
       - nan_percent = percentage of NaN residues in bfactor
       - chain_length = len(sequence)
  4. keep chains with nan_percent < 10, chain_length > 100 and bnorm_mean < 20
  5. drop (sequence, taxonomy_id) duplicates, keeping the best (lowest) resolution

Note: the source notebook computed nan_percent by grouping on pdb_id, which on an
array-valued column silently evaluates to 0 for every row (a no-op filter). We
apply the documented intent instead: the fraction of missing residues per chain.

Requires: pandas, numpy, pyarrow.

Examples:
    python preprocessing/clean_pdb.py \\
        --in data/processed/pdb_protein_features.parquet \\
        --out data/processed/pdb_protein_features_clean.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

XRAY_METHOD = "X-RAY DIFFRACTION"

# Chain-level quality cutoffs (see step 4).
MAX_NAN_PERCENT = 10.0  # at most 10% of residues may lack a B-factor
MIN_CHAIN_LENGTH = 100  # keep chains longer than this many residues
MAX_BNORM_MEAN = 20.0   # drop chains whose mean normalised B-factor is too high

# A (sequence, taxonomy_id) pair may appear in several entries; keep the sharpest.
IDENTITY_COLUMNS = ["sequence", "taxonomy_id"]

# Intermediate and metadata columns removed from the published output.
DROP_COLUMNS = [
    "method",
    "resolution",
    "wilson_b",
    "coverage",
    "bnorm_mean",
    "nan_percent",
]


def log(message: str) -> None:
    print(message, file=sys.stderr)


def nan_percent(bfactor: object) -> float:
    """Percentage of residues in a chain that have no B-factor."""
    values = np.asarray(bfactor, dtype=float)
    if values.size == 0:
        return 100.0
    return float(np.isnan(values).mean() * 100)


def clean(df: pd.DataFrame) -> pd.DataFrame:
    log(f"loaded: {len(df)} rows")

    # 1. Wilson B (and thus b_norm) only makes sense for X-ray structures.
    before = len(df)
    df = df[df["method"] == XRAY_METHOD]
    log(f"keep {XRAY_METHOD}: {before} -> {len(df)}")

    # 2. Need a positive Wilson B to normalise against.
    before = len(df)
    df = df[df["wilson_b"].notnull() & (df["wilson_b"] > 0)]
    log(f"keep usable wilson_b: {before} -> {len(df)}")

    df = df.copy()

    # 3. Normalise per-residue B-factors and derive chain-level statistics.
    df["bfactor_norm"] = df["bfactor"] / df["wilson_b"]
    df["bnorm_mean"] = df["bfactor_norm"].apply(np.nanmean)
    df["nan_percent"] = df["bfactor"].apply(nan_percent)
    df["length"] = df["sequence"].str.len()

    # 4. Keep only chains that are well measured, long enough and well ordered.
    before = len(df)
    df = df[
        (df["nan_percent"] < MAX_NAN_PERCENT)
        & (df["length"] > MIN_CHAIN_LENGTH)
        & (df["bnorm_mean"] < MAX_BNORM_MEAN)
    ]
    log(f"apply quality filters: {before} -> {len(df)}")

    # 5. Collapse duplicate sequences per organism to the best-resolved chain.
    before = len(df)
    df = df.sort_values("resolution", ascending=True).drop_duplicates(
        subset=IDENTITY_COLUMNS, keep="first"
    )
    log(f"drop duplicate {IDENTITY_COLUMNS}: {before} -> {len(df)}")

    # Drop intermediate/unneeded columns; keep only what the classifier consumes.
    df = df.drop(columns=DROP_COLUMNS)

    return df


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True, help="PDB features Parquet to clean")
    ap.add_argument("--out", required=True, help="cleaned Parquet to write")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    in_path = Path(args.inp)
    if not in_path.exists():
        log(f"missing input: {in_path}")
        return 1

    df = clean(pd.read_parquet(in_path))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, compression="zstd", index=False)

    log(f"done: {len(df)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
