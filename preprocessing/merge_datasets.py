#!/usr/bin/env python3
"""Merge the three processed protein datasets into one unified Parquet file.

The inputs describe proteins from different sources with different ids and
different per-residue features, so we union them on a common schema and keep a
`source` column to tell them apart:

    source         string        "disprot" | "mobidb" | "rcsb-pdb"
    id             string        source id: UniProt accession or PDB id
    organism       string        source organism; null if unknown
    taxonomy_id    string        NCBI taxonomy id; null if unknown
    sequence       string        one letter per residue
    disorder_mask  string        same length; '1' = disordered residue.
                                  Present for DisProt / MobiDB, null for PDB.
    coverage       string        same length; '1' = residue resolved in the
                                  structure. Present for PDB, null otherwise.
    bfactor        list<float32> same length; CA B-factor per residue.
                                  Present for PDB, null otherwise.

Each source contributes only the feature columns it has; the rest are null.
Data is streamed source by source so we never hold every dataset in memory.

Requires: pyarrow.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA = pa.schema(
    [
        ("source", pa.string()),
        ("id", pa.string()),
        # RCSB PDB only: the auth chain ids behind this row. With `id` (the PDB
        # entry) they form the key SIFTS maps to a UniProt accession; disprot and
        # mobidb carry an accession in `id` already, so theirs is null.
        ("chains", pa.list_(pa.string())),
        ("organism", pa.string()),
        ("taxonomy_id", pa.string()),
        ("sequence", pa.string()),
        ("disorder_mask", pa.string()),
        ("coverage", pa.string()),
        ("bfactor", pa.list_(pa.float32())),
    ]
)

# For each source: its label and how its id column maps onto our `id`.
# The feature columns (disorder_mask / coverage / bfactor) are taken by name
# when present and filled with nulls otherwise.
SOURCES = {
    "disprot": "Uniprot_ID",
    "mobidb": "Uniprot_ID",
    "rcsb-pdb": "pdb_id",
}


def normalize(table: pa.Table, source: str, id_column: str) -> pa.Table:
    """Reshape one source table onto SCHEMA, filling absent columns with nulls."""
    n = table.num_rows

    def column(name: str, dtype: pa.DataType) -> pa.Array:
        if name in table.column_names:
            return table.column(name).cast(dtype)
        return pa.nulls(n, dtype)

    return pa.table(
        {
            "source": pa.array([source] * n, pa.string()),
            "id": table.column(id_column).cast(pa.string()),
            "chains": column("chains", pa.list_(pa.string())),
            "organism": column("organism", pa.string()),
            "taxonomy_id": column("taxonomy_id", pa.string()),
            "sequence": column("sequence", pa.string()),
            "disorder_mask": column("disorder_mask", pa.string()),
            "coverage": column("coverage", pa.string()),
            "bfactor": column("bfactor", pa.list_(pa.float32())),
        },
        schema=SCHEMA,
    )


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--disprot", required=True, help="DisProt sequence/disorder Parquet")
    ap.add_argument("--mobidb", required=True, help="MobiDB Gold Parquet")
    ap.add_argument("--pdb", required=True, help="RCSB PDB protein features Parquet")
    ap.add_argument("--out", default="data/processed/merged_protein_dataset.parquet")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    inputs = {"disprot": args.disprot, "mobidb": args.mobidb, "rcsb-pdb": args.pdb}

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    with pq.ParquetWriter(out_path, SCHEMA, compression="zstd") as writer:
        for source, path in inputs.items():
            if not Path(path).exists():
                print(f"missing input for {source}: {path}", file=sys.stderr)
                return 1
            parquet = pq.ParquetFile(path)
            id_column = SOURCES[source]
            # stream row groups so a large source (MobiDB) never loads at once
            for group in parquet.iter_batches():
                table = normalize(pa.Table.from_batches([group]), source, id_column)
                writer.write_table(table)
            rows = parquet.metadata.num_rows
            total += rows
            print(f"  {source}: {rows} rows from {path}", file=sys.stderr)

    print(f"done: {total} rows -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
