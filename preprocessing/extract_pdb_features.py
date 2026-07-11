#!/usr/bin/env python3
"""Build a per-protein dataset from the local RCSB PDB mmCIF mirror.

Walks data/raw/pdb_mmCIF/**/<id>.cif.gz and writes one Parquet row per unique
protein (unique SEQRES sequence in a file). Columns:

    pdb_id         string        e.g. "101m"
    organism       string        source organism (gene source, not the
                                  expression host); null if unknown
    taxonomy_id    string        NCBI taxonomy id; null if unknown
    sequence       string        full SEQRES, one letter per residue
    disorder_mask  string        same length; '1' = residue in SEQRES but not
                                  resolved in the structure (disordered)
    bfactor        list<float32> same length; CA B-factor per residue,
                                  null where disorder_mask == '1'

Chains that share a sequence (homodimer copies) are merged into one row:
a residue counts as ordered if resolved in ANY copy, and its B-factor is the
min across copies. Chains with different sequences become separate rows.
Only the first model and only polypeptide(L) chains are used.

Requires: gemmi, pyarrow.
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import gemmi
import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA = pa.schema(
    [
        ("pdb_id", pa.string()),
        ("organism", pa.string()),
        ("taxonomy_id", pa.string()),
        ("sequence", pa.string()),
        ("disorder_mask", pa.string()),
        ("bfactor", pa.list_(pa.float32())),
    ]
)

# mmCIF source categories, in priority order: (category, organism_tag, taxid_tag).
# `_entity_src_gen` uses the *gene* source (natural origin), not the host.
SOURCE_CATEGORIES = [
    ("_entity_src_nat.", "pdbx_organism_scientific", "pdbx_ncbi_taxonomy_id"),
    ("_entity_src_gen.", "pdbx_gene_src_scientific_name", "pdbx_gene_src_ncbi_taxonomy_id"),
    ("_pdbx_entity_src_syn.", "organism_scientific", "ncbi_taxonomy_id"),
]


def read_source_organisms(block) -> dict[str, tuple[str, str | None]]:
    """Map entity_id -> (organism_name, taxonomy_id) from the mmCIF source blocks."""
    organisms: dict[str, tuple[str, str | None]] = {}
    for category, organism_tag, taxid_tag in SOURCE_CATEGORIES:
        try:
            table = block.find(category, ["entity_id", organism_tag, taxid_tag])
        except Exception:  # noqa: BLE001 - category or tag absent in this file
            continue
        for row in table:
            entity_id = row.str(0) if row.has(0) else None
            organism = row.str(1) if row.has(1) else None
            taxid = row.str(2) if row.has(2) else None
            if entity_id and organism and entity_id not in organisms:
                organisms[entity_id] = (organism, taxid)
    return organisms


def residue_bfactor(residue) -> float:
    """B-factor of the residue: the CA atom, or the min over atoms if there is no CA."""
    ca = residue.find_atom("CA", "*")
    if ca is not None:
        return round(ca.b_iso, 2)
    return round(min(atom.b_iso for atom in residue), 2)


def chain_residue_arrays(polymer, entity):
    """Return (sequence, mask, bfactor) aligned to the full SEQRES of one chain.

    mask[i] = 1 if residue i is missing from the coordinates (disordered),
    bfactor[i] = its CA B-factor, or None where it is missing. If the entity has
    no SEQRES, falls back to the modeled residues (all considered ordered).
    """
    seqres = entity.full_sequence if entity is not None else None
    if seqres:
        n = len(seqres)
        sequence = gemmi.one_letter_code(seqres)
        mask = [1] * n
        bfactor: list[float | None] = [None] * n
        for residue in polymer:
            i = residue.label_seq  # 1-based position in the SEQRES
            if i is None or not (1 <= i <= n):
                continue
            mask[i - 1] = 0
            bfactor[i - 1] = residue_bfactor(residue)
        return sequence, mask, bfactor

    sequence = polymer.make_one_letter_sequence()
    mask = [0] * len(sequence)
    bfactor = [residue_bfactor(residue) for residue in polymer]
    if len(bfactor) != len(sequence):
        return None
    return sequence, mask, bfactor


@dataclass
class Protein:
    """One unique protein sequence in a structure, with its per-residue arrays."""

    pdb_id: str
    organism: str | None
    taxonomy_id: str | None
    sequence: str
    mask: list[int]
    bfactor: list[float | None]

    def merge_copy(self, mask, bfactor) -> None:
        """Fold in another chain with the same sequence (a homodimer copy)."""
        for i, resolved in enumerate(mask):
            if resolved == 0:  # residue is present in this copy
                self.mask[i] = 0
                b = bfactor[i]
                if b is not None:
                    self.bfactor[i] = b if self.bfactor[i] is None else min(self.bfactor[i], b)


def extract_proteins(path: str):
    """Return (proteins, error) for one mmCIF file; error is a string or None."""
    pdb_id = Path(path).name.split(".")[0]
    try:
        block = gemmi.cif.read(path).sole_block()
        organisms = read_source_organisms(block)
        structure = gemmi.make_structure_from_block(block)
    except Exception as e:  # noqa: BLE001
        return [], f"{pdb_id}\t{type(e).__name__}: {e}"

    if len(structure) == 0:
        return [], None
    structure.setup_entities()

    by_sequence: dict[str, Protein] = {}
    for chain in structure[0]:  # first model only
        polymer = chain.get_polymer()
        if not polymer or polymer.check_polymer_type() != gemmi.PolymerType.PeptideL:
            continue
        entity = structure.get_entity_of(polymer)
        arrays = chain_residue_arrays(polymer, entity)
        if arrays is None:
            continue
        sequence, mask, bfactor = arrays
        if sequence in by_sequence:
            by_sequence[sequence].merge_copy(mask, bfactor)
        else:
            organism, taxid = organisms.get(entity.name, (None, None)) if entity else (None, None)
            by_sequence[sequence] = Protein(pdb_id, organism, taxid, sequence, mask, bfactor)

    return list(by_sequence.values()), None


def proteins_to_table(proteins: list[Protein]) -> pa.Table:
    """Build a Parquet table from a batch of proteins."""
    return pa.table(
        {
            "pdb_id": [p.pdb_id for p in proteins],
            "organism": [p.organism for p in proteins],
            "taxonomy_id": [p.taxonomy_id for p in proteins],
            "sequence": [p.sequence for p in proteins],
            "disorder_mask": ["".join(map(str, p.mask)) for p in proteins],
            "bfactor": pa.array([p.bfactor for p in proteins], type=pa.list_(pa.float32())),
        },
        schema=SCHEMA,
    )


def parse_args():
    ap = argparse.ArgumentParser(
        description="Per-protein dataset (sequence + disorder mask + CA B-factor) from RCSB PDB mmCIF."
    )
    ap.add_argument("--raw", default="data/raw/pdb_mmCIF", help="root of the RCSB PDB mmCIF mirror")
    ap.add_argument("--out", default="data/processed/pdb_protein_features.parquet")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="process only N files (debug)")
    ap.add_argument("--batch", type=int, default=20000, help="proteins per Parquet row group")
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    files = sorted(str(p) for p in Path(args.raw).rglob("*.cif.gz"))
    if args.limit:
        files = files[: args.limit]
    if not files:
        print(f"no .cif.gz under {args.raw}", file=sys.stderr)
        return 1
    print(f"processing {len(files)} files with {args.workers} workers ...", file=sys.stderr)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # buffer proteins and flush a Parquet row group once the batch is full,
    # so we never hold the whole dataset in memory
    batch: list[Protein] = []
    n_rows = n_err = done = 0

    with (
        pq.ParquetWriter(out_path, SCHEMA, compression="zstd") as writer,
        ProcessPoolExecutor(max_workers=args.workers) as pool,
    ):

        def flush():
            if batch:
                writer.write_table(proteins_to_table(batch))
                batch.clear()

        for proteins, error in pool.map(extract_proteins, files):
            done += 1
            if error:
                print(f"error: {error}", file=sys.stderr)
                n_err += 1
            batch.extend(proteins)
            n_rows += len(proteins)
            if len(batch) >= args.batch:
                flush()
            if done % 500 == 0 or done == len(files):
                print(f"  {done}/{len(files)} files | {n_rows} proteins | {n_err} errors",
                      file=sys.stderr)
        flush()

    print(f"done: {n_rows} proteins -> {args.out} ({n_err} errors)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
