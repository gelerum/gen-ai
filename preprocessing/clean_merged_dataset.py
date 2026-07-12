#!/usr/bin/env python3
"""Clean the merged protein dataset produced by merge_datasets.py.

The merge step unions DisProt, MobiDB and RCSB PDB row by row without touching
their values, so the merged Parquet still carries source-specific quirks: empty
strings that mean "missing", empty B-factor arrays, mixed-case organism names and
sequences, rows with no usable label, and exact duplicates contributed by more
than one source. This step normalises those and drops rows we cannot train on.

Steps (each logged with a before/after row count):
  1. empty strings -> null in every string column
  2. empty B-factor arrays -> null
  3. organism -> lower case, sequence -> upper case
  4. drop rows missing `sequence` or `taxonomy_id` (both are required keys)
  5. drop rows with no label at all (neither `disorder_mask` nor `bfactor`)
  6. drop exact duplicates on (taxonomy_id, sequence, disorder_mask, bfactor)

Whatever (sequence, taxonomy_id) duplicates remain after step 6 carry
*conflicting* labels across sources; we only count and report them here rather
than pick a winner, so the resolution policy stays an explicit later decision.

Requires: pandas, pyarrow.

Examples:
    python preprocessing/clean_merged_dataset.py \\
        --in data/processed/merged_protein_dataset.parquet \\
        --out data/processed/merged_protein_dataset_clean.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Columns that hold text; used for the empty-string -> null normalisation.
STRING_COLUMNS = [
    "source",
    "id",
    "organism",
    "taxonomy_id",
    "sequence",
    "disorder_mask",
    "coverage",
]

# A row is only useful if it has an identity (these keys) ...
REQUIRED_COLUMNS = ["sequence", "taxonomy_id"]
# ... and at least one per-residue label to learn from.
LABEL_COLUMNS = ["disorder_mask", "bfactor"]

# Exact-duplicate key: same organism + sequence + identical labels.
DEDUP_COLUMNS = ["taxonomy_id", "sequence", "disorder_mask", "bfactor"]
# Weaker key used only to report label conflicts that survive deduplication.
IDENTITY_COLUMNS = ["sequence", "taxonomy_id"]


def log(message: str) -> None:
    print(message, file=sys.stderr)


def is_empty_array(value: object) -> bool:
    return isinstance(value, (list, np.ndarray)) and len(value) == 0


def to_hashable(value: object) -> object:
    """Make a cell usable inside a duplicate key (arrays/lists -> tuples)."""
    if isinstance(value, np.ndarray):
        return tuple(value.tolist())
    if isinstance(value, list):
        return tuple(value)
    return value


def clean(df: pd.DataFrame) -> pd.DataFrame:
    log(f"loaded: {len(df)} rows")

    # 1. Empty strings are how some sources spell "missing"; make them null.
    present = [c for c in STRING_COLUMNS if c in df.columns]
    df[present] = df[present].replace("", pd.NA)

    # 2. An empty B-factor array carries no signal; treat it as missing.
    if "bfactor" in df.columns:
        df["bfactor"] = df["bfactor"].apply(lambda x: pd.NA if is_empty_array(x) else x)

    # 3. Canonicalise case so equal values compare equal during dedup.
    if "organism" in df.columns:
        df["organism"] = df["organism"].str.lower()
    if "sequence" in df.columns:
        df["sequence"] = df["sequence"].str.upper()

    # 4. Drop rows missing a required key.
    before = len(df)
    df = df.dropna(subset=REQUIRED_COLUMNS, how="any")
    log(f"drop missing {REQUIRED_COLUMNS}: {before} -> {len(df)}")

    # 5. Drop rows with no label at all.
    before = len(df)
    df = df.dropna(subset=LABEL_COLUMNS, how="all")
    log(f"drop rows with no label: {before} -> {len(df)}")

    # 6. Drop exact duplicates. Hash arrays into a throwaway key column so we
    #    can dedup on bfactor without mutating the column we write out.
    before = len(df)
    dedup_key = pd.DataFrame(
        {col: df[col].apply(to_hashable) for col in DEDUP_COLUMNS}
    )
    df = df[~dedup_key.duplicated(keep="first")]
    log(f"drop exact duplicates: {before} -> {len(df)}")

    # Report (do not resolve) remaining conflicting-label duplicates.
    conflicts = df.duplicated(subset=IDENTITY_COLUMNS, keep=False).sum()
    log(f"remaining {IDENTITY_COLUMNS} rows with conflicting labels: {conflicts}")

    return df


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True, help="merged Parquet to clean")
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
