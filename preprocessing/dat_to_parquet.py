#!/usr/bin/env python3
"""Parse a UniProt/SwissProt .dat.gz flat file into a Parquet table.

Reads the raw `uniprot_sprot.dat.gz` release downloaded by
pipeline/download_uniprot.nf and writes one row per entry, keeping the fields the
downstream steps need: `keywords`, the GO cross-references (`go`), and the other
database cross-references (Pfam, InterPro, AlphaFold, EC).

Entries are streamed and written in batches, so memory stays flat regardless of
release size.

Requires: biopython, pyarrow.

Examples:
    python preprocessing/dat_to_parquet.py \\
        --in data/raw/uniprot/uniprot_sprot.dat.gz \\
        --out data/processed/uniprot/uniprot_sprot.parquet
"""
from __future__ import annotations

import argparse
import gzip
import pyarrow as pa
import pyarrow.parquet as pq


from Bio import SwissProt
from pathlib import Path

BATCH_SIZE = 5000

def string_list(values):
    if values:
        return [str(x) for x in values]

    return []

def parse_gene_name(gene_name):
    if not gene_name:
        return None

    try:
        names = []
        for item in gene_name:

            if hasattr(item, "Name") and item.Name:
                names.append(item.Name)

        if names:
            return "; ".join(names)

    except Exception:
        pass

    return str(gene_name)

def parse_record(record):
    go, pdb, pfam, interpro, alphafold, ec_numbers = [], [], [], [], [], []

    for ref in record.cross_references:
        database = ref[0]

        if database == "GO":
            go.append(
                {
                    "id": str(ref[1]),
                    "term": str(ref[2]) if len(ref) > 2 else None,
                    "evidence": str(ref[3]) if len(ref) > 3 else None,
                }
            )

        elif database == "PDB":
            pdb.append(str(ref[1]))

        elif database == "Pfam":
            pfam.append(str(ref[1]))

        elif database == "InterPro":
            interpro.append(str(ref[1]))

        elif database == "AlphaFoldDB":
            alphafold.append(str(ref[1]))

        elif database == "EC":
            ec_numbers.append(str(ref[1]))

    return {
        "accession":
            str(record.accessions[0]),
        "accessions":
            string_list(record.accessions),
        "entry_name":
            str(record.entry_name),
        "reviewed":
            True,
        "protein_name":
            str(record.description),
        "gene_name":
            parse_gene_name(record.gene_name),
        "organism":
            " ".join(record.organism),
        "taxonomy":
            string_list(record.organism_classification),
        "taxonomy_id":
            string_list(record.taxonomy_id),
        "host_organism":
            string_list(record.host_organism),
        "keywords":
            string_list(record.keywords),
        "sequence":
            str(record.sequence),
        "length":
            int(len(record.sequence)),
        "references":
            int(len(record.references)),
        "go":
            go,
        "pdb":
            pdb,
        "pfam":
            pfam,
        "interpro":
            interpro,
        "alphafold":
            alphafold,
        "ec_numbers":
            ec_numbers,
    }

GO_TYPE = pa.struct(
    [
        ("id", pa.string()),
        ("term", pa.string()),
        ("evidence", pa.string()),
    ]
)


SCHEMA = pa.schema(
    [
        ("accession", pa.string()),
        ("accessions", pa.list_(pa.string())),
        ("entry_name", pa.string()),
        ("reviewed", pa.bool_()),
        ("protein_name", pa.string()),
        ("gene_name", pa.string()),
        ("organism", pa.string()),
        ("taxonomy", pa.list_(pa.string())),
        ("taxonomy_id", pa.list_(pa.string())),
        ("host_organism", pa.list_(pa.string())),
        ("keywords", pa.list_(pa.string())),
        ("sequence", pa.string()),
        ("length", pa.int64()),
        ("references", pa.int64()),
        ("go", pa.list_(GO_TYPE)),
        ("pdb", pa.list_(pa.string())),
        ("pfam", pa.list_(pa.string())),
        ("interpro", pa.list_(pa.string())),
        ("alphafold", pa.list_(pa.string())),
        ("ec_numbers", pa.list_(pa.string())),
    ]
)

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="input", required=True, help="uniprot_sprot .dat.gz to read")
    ap.add_argument("--out", required=True, help="output .parquet path")

    return ap.parse_args()

def main():
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.out)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    writer = pq.ParquetWriter(
        output_path,
        SCHEMA,
        compression="zstd",
    )

    batch = []
    counter = 0

    with gzip.open(input_path, "rt") as handle:
        for record in SwissProt.parse(handle):
            batch.append(parse_record(record))

            if len(batch) >= BATCH_SIZE:
                table = pa.Table.from_pylist(batch, schema=SCHEMA)

                writer.write_table(table)

                counter += len(batch)

                print(
                    f"{counter:,} records parsed"
                )

                batch.clear()

    if batch:
        table = pa.Table.from_pylist(batch, schema=SCHEMA)
        writer.write_table(table)
        counter += len(batch)

    writer.close()

if __name__ == "__main__":
    main()