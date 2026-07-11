#!/usr/bin/env python3
"""Build the DisProt sequence/disorder-mask dataset.

The Nextflow pipeline handles orchestration and downloads. This script owns the
data transformations:

* extract one row per UniProt accession from the DisProt TSV;
* build a sequence-aligned 0/1 disorder mask for one FASTA file;
* merge row shards and write the final Parquet dataset plus error report.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


DATASET_SCHEMA = ["Uniprot_ID", "organism", "taxonomy_id", "sequence", "disorder_mask"]
ERROR_SCHEMA = ["Uniprot_ID", "error_type", "detail"]


def parse_regions(value: str) -> list[int]:
    """Parse a comma-separated list of 1-based region coordinates."""
    return [int(item) for item in value.split(",") if item and item != "null"]


def read_fasta_sequence(path: Path) -> str:
    """Read a FASTA file and return the concatenated sequence."""
    sequence_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                sequence_parts.append(line)
    return "".join(sequence_parts)


def extract_proteins(args: argparse.Namespace) -> int:
    """Extract unique UniProt proteins and disorder intervals from DisProt TSV."""
    proteins: dict[str, dict[str, object]] = {}

    with Path(args.disprot_tsv).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            acc = row["UniProt ACC"].strip()
            if not acc:
                continue

            item = proteins.setdefault(
                acc,
                {
                    "Uniprot_ID": acc,
                    "organism": row["Organism"],
                    "taxonomy_id": row["NCBI Taxon ID"],
                    "starts": [],
                    "ends": [],
                },
            )

            is_disorder = row["Term namespace"] == "Structural state" and row["Term name"] == "disorder"
            if not is_disorder:
                continue

            try:
                start = int(row["Start"])
                end = int(row["End"])
            except ValueError:
                continue

            if start > 0 and end >= start:
                item["starts"].append(start)
                item["ends"].append(end)

    rows = sorted(proteins.values(), key=lambda item: item["Uniprot_ID"])
    if args.limit > 0:
        rows = rows[: args.limit]

    with Path(args.out).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Uniprot_ID", "organism", "taxonomy_id", "starts", "ends"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "Uniprot_ID": row["Uniprot_ID"],
                    "organism": row["organism"],
                    "taxonomy_id": row["taxonomy_id"],
                    "starts": ",".join(map(str, row["starts"])),
                    "ends": ",".join(map(str, row["ends"])),
                }
            )

    return 0


def build_row(args: argparse.Namespace) -> int:
    """Build one dataset row, or one error row, for a UniProt FASTA file."""
    row_out = Path(args.row_out)
    error_out = Path(args.error_out)
    row_out.write_text("")
    error_out.write_text("")

    sequence = read_fasta_sequence(Path(args.fasta))
    if not sequence:
        with error_out.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([args.acc, "empty_fasta", args.fasta])
        return 0

    starts = parse_regions(args.starts)
    ends = parse_regions(args.ends)
    mask = ["0"] * len(sequence)

    for start, end in zip(starts, ends):
        left = max(start, 1) - 1
        right = min(end, len(sequence))
        for idx in range(left, right):
            mask[idx] = "1"

    with row_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([args.acc, args.organism, args.taxonomy_id, sequence, "".join(mask)])

    return 0


def write_parquet(args: argparse.Namespace) -> int:
    """Merge row shards and write Parquet plus a TSV error report."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        print(
            "pyarrow is required to write Parquet. Install it in the Nextflow environment, "
            "for example: conda install -c conda-forge pyarrow",
            file=sys.stderr,
        )
        raise

    rows = []
    for row_file in sorted(Path(args.rows_dir).glob("*.dataset_row.tsv")):
        with row_file.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for values in reader:
                if not values:
                    continue
                rows.append(
                    {
                        "Uniprot_ID": values[0],
                        "organism": values[1],
                        "taxonomy_id": values[2],
                        "sequence": values[3],
                        "disorder_mask": values[4],
                    }
                )

    rows.sort(key=lambda item: item["Uniprot_ID"])

    table = pa.Table.from_pylist(
        rows,
        schema=pa.schema(
            [
                ("Uniprot_ID", pa.string()),
                ("organism", pa.string()),
                ("taxonomy_id", pa.string()),
                ("sequence", pa.string()),
                ("disorder_mask", pa.string()),
            ]
        ),
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path)

    errors = []
    for error_file in sorted(Path(args.errors_dir).glob("*.dataset_error.tsv")):
        with error_file.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for values in reader:
                if not values:
                    continue
                errors.append(
                    {
                        "Uniprot_ID": values[0],
                        "error_type": values[1],
                        "detail": values[2] if len(values) > 2 else "",
                    }
                )

    error_path = Path(args.errors_out)
    error_path.parent.mkdir(parents=True, exist_ok=True)
    with error_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ERROR_SCHEMA, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(errors)

    with Path(args.manifest).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "file", "rows"])
        writer.writerow(["disprot_parquet", str(out_path), len(rows)])
        writer.writerow(["disprot_errors", str(error_path), len(errors)])

    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract-proteins")
    extract.add_argument("--disprot-tsv", required=True)
    extract.add_argument("--out", required=True)
    extract.add_argument("--limit", type=int, default=0)
    extract.set_defaults(func=extract_proteins)

    row = subparsers.add_parser("build-row")
    row.add_argument("--acc", required=True)
    row.add_argument("--organism", required=True)
    row.add_argument("--taxonomy-id", required=True)
    row.add_argument("--starts", default="")
    row.add_argument("--ends", default="")
    row.add_argument("--fasta", required=True)
    row.add_argument("--row-out", required=True)
    row.add_argument("--error-out", required=True)
    row.set_defaults(func=build_row)

    parquet = subparsers.add_parser("write-parquet")
    parquet.add_argument("--rows-dir", default=".")
    parquet.add_argument("--errors-dir", default=".")
    parquet.add_argument("--out", required=True)
    parquet.add_argument("--errors-out", required=True)
    parquet.add_argument("--manifest", default="disprot_dataset.manifest.tsv")
    parquet.set_defaults(func=write_parquet)

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
