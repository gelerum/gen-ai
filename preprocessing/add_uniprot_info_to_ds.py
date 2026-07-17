#!/usr/bin/env python3
"""Attach UniProt `keywords` and `go` to the merged protein dataset.

The merged dataset unions three sources that identify proteins differently, so
each needs its own route to a UniProt accession:

  * disprot / mobidb -- `id` already *is* a UniProt accession; joined directly.
  * rcsb-pdb         -- `id` is a PDB entry id (e.g. `101m`); resolved through
                        the SIFTS PDB<->UniProt mapping.

A PDB entry generally holds several chains that map to several *different*
accessions, so joining on the entry id alone is a cartesian product: it inflates
the rcsb-pdb rows ~13x and silently duplicates labels. There are two routes out,
picked by whether the dataset carries the `chains` column:

  * `chains` present -- exact. SIFTS keys on (PDB, CHAIN) and `chains` holds the
    auth chain ids behind the row, so each row resolves to its own accession.
  * `chains` absent  -- fallback. Annotations are aggregated per PDB entry, so
    every row of an entry gets the union over all of its accessions. A
    heteromeric complex (4hhb: alpha + beta) then smears the annotations of its
    members across both rows.

Either way the row count is preserved exactly, and the script refuses to write
if it is not. Regenerate pdb_protein_features.parquet with a current
extract_pdb_features.py to get `chains` and the exact route.

Requires: pandas, pyarrow.

Examples:
    python preprocessing/add_uniprot_info_to_ds.py \\
        --in data/processed/merged_protein_dataset.parquet \\
        --uniprot data/processed/uniprot/uniprot_sprot.parquet \\
        --sifts data/raw/sifts/pdb_chain_uniprot.csv.gz \\
        --out data/processed/merged_dataset_with_uniprot.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Sources whose `id` is already a UniProt accession.
ACCESSION_SOURCES = ["disprot", "mobidb"]
# Source whose `id` is a PDB entry that needs the SIFTS mapping.
PDB_SOURCE = "rcsb-pdb"

ANNOTATION_COLUMNS = ["keywords", "go"]


def log(message):
    print(message, file=sys.stderr, flush=True)


def iter_lists(series):
    """Yield only the real list values of a group.

    An accession that SIFTS knows but SwissProt does not (a TrEMBL entry) comes
    out of the left join as NaN rather than None, so test for the list itself.
    """
    for values in series:
        if isinstance(values, (list, np.ndarray)):
            yield values


def unique_keywords(series):
    """Flatten a group of keyword lists into one order-preserving unique list."""
    seen = {}

    for values in iter_lists(series):
        for value in values:
            seen[value] = None

    return list(seen)


def unique_go(series):
    """Flatten a group of GO lists, keeping the first record per GO id."""
    seen = {}

    for values in iter_lists(series):
        for value in values:
            go_id = value["id"]

            if go_id not in seen:
                seen[go_id] = value

    return list(seen.values())


def aggregate_by_pdb(mapping, uniprot):
    """Build one keywords/go record per PDB entry from its mapped accessions."""
    annotated = mapping.merge(uniprot, on="accession", how="left")

    return annotated.groupby("pdb_id").agg(
        keywords=("keywords", unique_keywords),
        go=("go", unique_go),
    )


def annotate_by_chain(by_pdb, mapping, uniprot):
    """Exact route: resolve each row's own chains to one accession via SIFTS.

    `chains` lists the auth chain ids folded into the row, which is exactly what
    SIFTS keys on. All chains of a row hold the same sequence and so resolve to
    the same accession; the first non-null one is taken.
    """
    exploded = by_pdb[["id", "chains"]].explode("chains")
    exploded = exploded.merge(
        mapping,
        left_on=["id", "chains"],
        right_on=["pdb_id", "chain"],
        how="left",
    ).set_index(exploded.index)

    # groupby.first() skips nulls, so a row keeps its accession even when some of
    # its chains are absent from SIFTS.
    accession = exploded.groupby(level=0)["accession"].first()

    annotated = by_pdb.assign(accession=accession).merge(uniprot, on="accession", how="left")

    return annotated.drop(columns="accession")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="input", required=True, help="merged dataset Parquet to annotate")
    ap.add_argument("--uniprot", required=True, help="uniprot_sprot Parquet from dat_to_parquet.py")
    ap.add_argument("--sifts", required=True, help="SIFTS pdb_chain_uniprot.csv.gz")
    ap.add_argument("--out", required=True, help="annotated Parquet to write")

    return ap.parse_args()


def main():
    args = parse_args()

    output_path = Path(args.out)

    df = pd.read_parquet(args.input)
    log(f"dataset: {len(df):,} rows")

    uniprot = pd.read_parquet(args.uniprot, columns=["accession"] + ANNOTATION_COLUMNS)
    uniprot = uniprot.drop_duplicates(subset="accession")
    log(f"uniprot: {len(uniprot):,} accessions")

    sifts = pd.read_csv(
        args.sifts,
        compression="gzip",
        comment="#",
        usecols=["PDB", "CHAIN", "SP_PRIMARY"],
        low_memory=False,
    )
    sifts = sifts.drop_duplicates().rename(
        columns={"PDB": "pdb_id", "CHAIN": "chain", "SP_PRIMARY": "accession"}
    )
    log(f"sifts: {len(sifts):,} rows over {sifts['pdb_id'].nunique():,} PDB entries")

    unknown = df[~df["source"].isin(ACCESSION_SOURCES + [PDB_SOURCE])]
    if len(unknown):
        raise SystemExit(f"no accession route for sources: {sorted(unknown['source'].unique())}")

    by_accession = df[df["source"].isin(ACCESSION_SOURCES)]
    by_pdb = df[df["source"] == PDB_SOURCE]

    # disprot / mobidb: `id` is the accession, so a plain left join stays 1:1.
    annotated_by_accession = by_accession.merge(
        uniprot.rename(columns={"accession": "id"}),
        on="id",
        how="left",
    )

    # rcsb-pdb: use the exact per-chain route when the dataset carries `chains`,
    # otherwise fall back to smearing each entry's annotations over its rows.
    if "chains" in by_pdb.columns and by_pdb["chains"].notna().any():
        log("rcsb-pdb: exact (PDB, CHAIN) join via SIFTS")
        annotated_by_pdb = annotate_by_chain(by_pdb, sifts, uniprot)
    else:
        log("rcsb-pdb: no `chains` column, falling back to per-PDB aggregation")
        mapping = sifts[["pdb_id", "accession"]].drop_duplicates()
        aggregated = aggregate_by_pdb(mapping, uniprot)
        annotated_by_pdb = by_pdb.merge(aggregated, left_on="id", right_index=True, how="left")

    merged = pd.concat([annotated_by_accession, annotated_by_pdb], ignore_index=True)

    if len(merged) != len(df):
        raise SystemExit(f"row count changed: {len(df):,} -> {len(merged):,}, refusing to write")

    # A PDB entry that SIFTS maps only to non-SwissProt accessions aggregates to
    # an empty list rather than null, so count by length, not by notna().
    annotated = merged["keywords"].apply(lambda v: isinstance(v, (list, np.ndarray)) and len(v) > 0)

    for source in sorted(merged["source"].unique()):
        rows = merged["source"] == source
        hits = int(annotated[rows].sum())
        log(f"{source:10s} {int(rows.sum()):>9,} rows | {hits:>9,} annotated ({hits / int(rows.sum()):.1%})")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(output_path, compression="zstd", index=False)
    log(f"wrote {output_path} ({len(merged):,} rows)")


if __name__ == "__main__":
    main()
