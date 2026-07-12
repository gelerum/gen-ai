#!/usr/bin/env python3
"""Convert a MobiDB Gold .mjson.gz archive to Parquet.

The source archive is gzip-compressed JSON Lines: each line is one MobiDB
protein record. This script keeps processing streaming-friendly and can write
either:

* a compact sequence/disorder dataset with one typed row per protein;
* a full dump keeping each original record as a JSON string column.

Output is a single Parquet file written in batches so the whole archive never
has to be held in memory.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import sys
from pathlib import Path
from typing import Any, Iterator

import pyarrow as pa
import pyarrow.parquet as pq


DEFAULT_DISORDER_KEY = "curated-disorder-merge"
DEFAULT_PRIMARY_VARIANT = "curated"
DEFAULT_DISORDER_VARIANTS = (
    "curated=curated-disorder-merge;"
    "homology=homology-disorder-merge;"
    "prediction_priority=prediction-disorder-priority;"
    "prediction_mobidb_lite=prediction-disorder-mobidb_lite;"
    "all_priority=curated-disorder-merge,homology-disorder-merge,prediction-disorder-priority"
)


def log_progress(message: str, progress_log: Path | None = None) -> None:
    """Print progress and optionally append it to a persistent log file."""
    line = f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {message}"
    print(line, flush=True)
    if progress_log is not None:
        progress_log.parent.mkdir(parents=True, exist_ok=True)
        with progress_log.open("a", encoding="utf-8") as handle:
            handle.write(f"{line}\n")


def open_maybe_gzip(path: Path):
    """Open plain text or gzip-compressed input for line-wise reading."""
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def iter_records(path: Path) -> Iterator[dict[str, Any]]:
    """Yield JSON objects from a MobiDB JSON Lines archive."""
    with open_maybe_gzip(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if not isinstance(record, dict):
                raise ValueError(f"{path}:{line_number}: expected a JSON object")
            yield record


def mongo_oid(value: Any) -> str | None:
    """Return the Mongo ObjectId string from a MobiDB _id field when present."""
    if isinstance(value, dict) and "$oid" in value:
        return str(value["$oid"])
    if value is None:
        return None
    return str(value)


def normalize_taxonomy(value: Any) -> str:
    """Normalize MobiDB taxonomy fields to a stable string."""
    if value is None:
        return ""
    if isinstance(value, dict):
        for key in ("id", "taxid", "ncbi_taxon_id", "ncbi_taxonomy_id"):
            if key in value and value[key] is not None:
                return str(value[key])
        return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return str(value)


def normalize_regions(value: Any) -> list[list[int]]:
    """Return valid 1-based inclusive regions from an annotation block."""
    if not isinstance(value, dict):
        return []

    regions = value.get("regions", [])
    if not isinstance(regions, list):
        return []

    normalized: list[list[int]] = []
    for region in regions:
        if not isinstance(region, (list, tuple)) or len(region) < 2:
            continue
        try:
            start = int(region[0])
            end = int(region[1])
        except (TypeError, ValueError):
            continue
        if start > 0 and end >= start:
            normalized.append([start, end])

    normalized.sort(key=lambda item: (item[0], item[1]))
    return normalized


def infer_length(record: dict[str, Any]) -> int:
    """Infer protein length from explicit length, sequence, or score arrays."""
    try:
        length = int(record.get("length") or 0)
    except (TypeError, ValueError):
        length = 0

    sequence = record.get("sequence")
    if isinstance(sequence, str) and len(sequence) > length:
        length = len(sequence)

    if length > 0:
        return length

    for value in record.values():
        if isinstance(value, dict) and isinstance(value.get("scores"), list):
            return len(value["scores"])

    return 0


def build_mask(length: int, regions: list[list[int]]) -> str:
    """Build a 0/1 mask from 1-based inclusive regions."""
    if length <= 0:
        return ""

    mask = ["0"] * length
    for start, end in regions:
        left = max(start, 1) - 1
        right = min(end, length)
        for idx in range(left, right):
            mask[idx] = "1"
    return "".join(mask)


def merge_regions(regions: list[list[int]]) -> list[list[int]]:
    """Merge overlapping or adjacent 1-based inclusive regions."""
    if not regions:
        return []

    ordered = sorted(regions, key=lambda item: (item[0], item[1]))
    merged = [ordered[0][:]]
    for start, end in ordered[1:]:
        last = merged[-1]
        if start <= last[1] + 1:
            last[1] = max(last[1], end)
        else:
            merged.append([start, end])
    return merged


def parse_disorder_keys(value: str) -> list[str]:
    """Parse comma-separated MobiDB disorder annotation keys."""
    keys = [item.strip() for item in value.split(",") if item.strip()]
    return keys or [DEFAULT_DISORDER_KEY]


def parse_disorder_variants(value: str) -> dict[str, list[str]]:
    """Parse variant specs like name=key1,key2;other=key3."""
    variants: dict[str, list[str]] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            name, keys_value = item.split("=", 1)
            name = name.strip()
        else:
            keys_value = item
            name = item.replace("-", "_").replace(",", "_")
        keys = parse_disorder_keys(keys_value)
        if name:
            variants[name] = keys

    return variants or parse_disorder_variants(DEFAULT_DISORDER_VARIANTS)


def collect_disorder_regions(record: dict[str, Any], disorder_keys: list[str]) -> tuple[list[list[int]], list[str]]:
    """Collect disorder regions from one or more annotation blocks."""
    collected: list[list[int]] = []
    sources: list[str] = []

    for key in disorder_keys:
        regions = normalize_regions(record.get(key))
        if not regions:
            continue

        collected.extend(regions)
        sources.append(key)

    return merge_regions(collected), sources


def build_disorder_variant(record: dict[str, Any], length: int, disorder_keys: list[str]) -> dict[str, Any]:
    """Build one named disorder-mask variant."""
    regions, sources = collect_disorder_regions(record, disorder_keys)
    mask = build_mask(length, regions)
    content_count = mask.count("1")
    return {
        "sources": sources,
        "regions": regions,
        "content_count": content_count,
        "content_fraction": round(content_count / length, 6) if length else 0,
        "mask": mask,
    }


def compact_record(
    record: dict[str, Any],
    disorder_variants: dict[str, list[str]],
    primary_disorder_variant: str,
) -> dict[str, Any]:
    """Build the compact MobiDB dataset row."""
    sequence = record.get("sequence")
    if not isinstance(sequence, str):
        sequence = ""

    length = infer_length(record)
    disorder_masks = {
        name: build_disorder_variant(record, length, keys)
        for name, keys in disorder_variants.items()
    }

    primary_name = primary_disorder_variant
    if primary_name not in disorder_masks:
        primary_name = next(iter(disorder_masks))
    primary = disorder_masks[primary_name]

    return {
        "Uniprot_ID": str(record.get("acc") or ""),
        "mobidb_id": mongo_oid(record.get("_id")),
        "organism": str(record.get("organism") or ""),
        "taxonomy_id": normalize_taxonomy(record.get("ncbi_taxon_id") or record.get("taxonomy")),
        "length": length,
        "sequence": sequence,
        "primary_disorder_variant": primary_name,
        "disorder_masks": disorder_masks,
        "disorder_source": ",".join(disorder_variants[primary_name]),
        "disorder_sources": primary["sources"],
        "disorder_regions": primary["regions"],
        "disorder_content_count": primary["content_count"],
        "disorder_content_fraction": primary["content_fraction"],
        "disorder_mask": primary["mask"],
    }


def dataset_schema(variant_names: list[str]) -> pa.Schema:
    """Build the Arrow schema for the compact per-protein dataset."""
    variant_struct = pa.struct(
        [
            ("sources", pa.list_(pa.string())),
            ("regions", pa.list_(pa.list_(pa.int32()))),
            ("content_count", pa.int32()),
            ("content_fraction", pa.float64()),
            ("mask", pa.string()),
        ]
    )
    masks_struct = pa.struct([(name, variant_struct) for name in variant_names])
    return pa.schema(
        [
            ("Uniprot_ID", pa.string()),
            ("mobidb_id", pa.string()),
            ("organism", pa.string()),
            ("taxonomy_id", pa.string()),
            ("length", pa.int32()),
            ("sequence", pa.string()),
            ("primary_disorder_variant", pa.string()),
            ("disorder_masks", masks_struct),
            ("disorder_source", pa.string()),
            ("disorder_sources", pa.list_(pa.string())),
            ("disorder_regions", pa.list_(pa.list_(pa.int32()))),
            ("disorder_content_count", pa.int32()),
            ("disorder_content_fraction", pa.float64()),
            ("disorder_mask", pa.string()),
        ]
    )


def full_schema() -> pa.Schema:
    """Schema for the full dump: one JSON string per original record."""
    return pa.schema([("record", pa.string())])


def write_parquet(
    records: Iterator[dict[str, Any]],
    out_path: Path,
    schema: pa.Schema,
    progress_log: Path | None = None,
    progress_every: int = 10_000,
    batch_size: int = 2_000,
) -> int:
    """Stream records into a single Parquet file, one batch at a time."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    batch: list[dict[str, Any]] = []
    writer = pq.ParquetWriter(out_path, schema, compression="zstd")
    try:
        for record in records:
            batch.append(record)
            count += 1
            if len(batch) >= batch_size:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                batch.clear()
            if progress_every > 0 and count % progress_every == 0:
                log_progress(f"[convertMobidbGoldToJson] converted {count} records", progress_log)
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
    finally:
        writer.close()
    log_progress(f"[convertMobidbGoldToJson] converted {count} records total", progress_log)
    return count


def write_manifest(path: Path, source: Path, out: Path, mode: str, rows: int, disorder_variants: dict[str, list[str]]) -> None:
    """Write a small TSV manifest next to the Nextflow task output."""
    with path.open("w", encoding="utf-8") as handle:
        handle.write("source\tfile\tmode\tdisorder_variants\trows\n")
        variants_json = json.dumps(disorder_variants, sort_keys=True, separators=(",", ":"))
        handle.write(f"mobidb_gold\t{out}\t{mode}\t{variants_json}\t{rows}\n")
        handle.write(f"mobidb_raw\t{source}\traw\t\t\n")


def convert(args: argparse.Namespace) -> int:
    in_path = Path(args.input)
    out_path = Path(args.out)
    manifest_path = Path(args.manifest)
    progress_log = Path(args.progress_log) if args.progress_log else None
    variants_value = args.disorder_variants
    if not variants_value and (args.disorder_keys or args.disorder_key != DEFAULT_DISORDER_KEY):
        variants_value = f"{DEFAULT_PRIMARY_VARIANT}={args.disorder_keys or args.disorder_key}"
    disorder_variants = parse_disorder_variants(variants_value or DEFAULT_DISORDER_VARIANTS)

    if args.mode == "full":
        schema = full_schema()
        records = (
            {"record": json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":"))}
            for record in iter_records(in_path)
        )
    else:
        schema = dataset_schema(list(disorder_variants.keys()))
        records = (
            compact_record(record, disorder_variants, args.primary_disorder_variant)
            for record in iter_records(in_path)
        )

    log_progress(f"[convertMobidbGoldToJson] converting {in_path} to {out_path}", progress_log)
    rows = write_parquet(records, out_path, schema, progress_log=progress_log)
    write_manifest(manifest_path, in_path, out_path, args.mode, rows, disorder_variants)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="MobiDB .mjson or .mjson.gz archive")
    parser.add_argument("--out", required=True, help="output .parquet path")
    parser.add_argument("--manifest", default="mobidb_gold_json.manifest.tsv")
    parser.add_argument(
        "--mode",
        choices=("dataset", "full"),
        default="dataset",
        help="dataset writes compact protein rows; full writes original records as a JSON array",
    )
    parser.add_argument(
        "--disorder-key",
        default=DEFAULT_DISORDER_KEY,
        help="single MobiDB annotation block used for the compact disorder mask",
    )
    parser.add_argument(
        "--disorder-keys",
        default="",
        help="comma-separated MobiDB annotation blocks used for the compact disorder mask",
    )
    parser.add_argument(
        "--disorder-variants",
        default="",
        help="semicolon-separated variants, for example curated=curated-disorder-merge;all=key1,key2",
    )
    parser.add_argument(
        "--primary-disorder-variant",
        default=DEFAULT_PRIMARY_VARIANT,
        help="variant copied to legacy top-level disorder_mask fields",
    )
    parser.add_argument("--progress-log", default="", help="append conversion progress to this log file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        return convert(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
