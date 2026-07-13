#!/usr/bin/env python3
"""Attach ESM2 embeddings to any Parquet with a `sequence` column.

General-purpose embedder: point it at a Parquet file and it adds one fixed-size
`sequence_embedding` vector per row, produced by an ESM2 protein language model
(``facebook/esm2_t6_8M_UR50D`` by default, hidden size 320). It does not filter
or reshape the input — every row is kept, only the new column is added.

The embedding is a mean pool over the model's per-residue hidden states, with the
special ``<cls>``/``<eos>`` and padding tokens masked out, so each sequence maps
to a single vector regardless of length. Long sequences are split into
non-overlapping windows and mean-pooled across all of them, so nothing is
truncated. Identical sequences are embedded once and the vector is broadcast to
every row that carries it.

The model is downloaded from the Hugging Face Hub on first use; it runs on CUDA
when available and falls back to CPU otherwise.

Requires: pandas, pyarrow, numpy, torch, transformers.

Examples:
    python preprocessing/embed_sequences.py \\
        --in data/processed/pdb_protein_features_clean.parquet \\
        --out data/processed/pdb_protein_features_clean_esm2.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Column that holds the amino-acid string, and the column we write.
SEQUENCE_COLUMN = "sequence"
EMBEDDING_COLUMN = "sequence_embedding"

# Self-attention memory grows as O(L^2), so we cap each forward pass at this many
# residues; longer sequences are split into windows and mean-pooled across them.
WINDOW_RESIDUES = 1022


def log(message: str) -> None:
    print(message, file=sys.stderr)


def windows(seq: str, window: int) -> list[str]:
    """Split a sequence into consecutive non-overlapping residue windows."""
    if len(seq) <= window:
        return [seq]
    return [seq[i : i + window] for i in range(0, len(seq), window)]


def embed_sequences(
    sequences: list[str],
    model_name: str,
    batch_size: int,
    window: int,
    device: str | None,
) -> np.ndarray:
    """Return an (N, hidden) float32 array, one mean-pooled vector per sequence.

    Each sequence is split into windows of at most `window` residues; every
    window is embedded and reduced to a residue sum + count, and the per-window
    sums/counts are accumulated so the final vector is the mean over *all*
    residues of the sequence. Windows are processed longest-first so each padded
    batch wastes as little compute as possible.
    """
    # Imported here so the module stays importable without the heavy ML stack
    # and so nothing is downloaded until this stage is actually executed.
    import torch
    from transformers import AutoModel, AutoTokenizer

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"loading {model_name} on {device}")

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval().to(device)

    # Half precision is a free speed/memory win on GPU; keep fp32 on CPU.
    autocast_dtype = torch.float16 if device == "cuda" else None

    hidden_size = model.config.hidden_size
    # Accumulate residue sums and counts per sequence across its windows.
    sums = np.zeros((len(sequences), hidden_size), dtype=np.float64)
    counts = np.zeros(len(sequences), dtype=np.int64)

    # One work item per window, tagged with the sequence it belongs to.
    items = [
        (seq_i, chunk)
        for seq_i, seq in enumerate(sequences)
        for chunk in windows(seq, window)
    ]
    n_long = sum(1 for seq in sequences if len(seq) > window)
    log(f"{len(items)} windows from {len(sequences)} sequences ({n_long} split across >1 window)")

    # Longest-first ordering minimises padding within each batch.
    items.sort(key=lambda it: len(it[1]), reverse=True)

    special_ids = {
        tid
        for tid in (
            tokenizer.cls_token_id,
            tokenizer.eos_token_id,
            tokenizer.pad_token_id,
            tokenizer.bos_token_id,
            tokenizer.sep_token_id,
        )
        if tid is not None
    }

    total = len(items)
    with torch.no_grad():
        for start in range(0, total, batch_size):
            batch = items[start : start + batch_size]
            batch_seqs = [chunk for _, chunk in batch]

            enc = tokenizer(
                batch_seqs,
                padding=True,
                return_tensors="pt",
            ).to(device)

            if autocast_dtype is not None:
                with torch.autocast(device_type="cuda", dtype=autocast_dtype):
                    hidden = model(**enc).last_hidden_state
            else:
                hidden = model(**enc).last_hidden_state

            # Keep only real residue tokens: attended AND not a special token.
            keep = enc["attention_mask"].clone().bool()
            for tid in special_ids:
                keep &= enc["input_ids"] != tid
            keep_f = keep.unsqueeze(-1)  # (B, T, 1)

            win_sum = (hidden * keep_f).sum(dim=1).float().cpu().numpy()
            win_count = keep.sum(dim=1).cpu().numpy()

            for local, (seq_i, _) in enumerate(batch):
                sums[seq_i] += win_sum[local]
                counts[seq_i] += win_count[local]

            done = min(start + batch_size, total)
            log(f"embedded {done}/{total} windows")

    return (sums / np.clip(counts, 1, None)[:, None]).astype(np.float32)


def run(
    in_path: Path,
    model_name: str,
    batch_size: int,
    window: int,
    device: str | None,
) -> pd.DataFrame:
    df = pd.read_parquet(in_path)
    log(f"loaded: {len(df)} rows")
    if SEQUENCE_COLUMN not in df.columns:
        raise SystemExit(f"input has no '{SEQUENCE_COLUMN}' column")

    # The embedding depends only on the sequence, so run the model once per
    # unique sequence and broadcast the vectors back onto every row.
    unique_seqs = df[SEQUENCE_COLUMN].dropna().drop_duplicates().tolist()
    log(f"unique sequences: {len(unique_seqs)} (of {len(df)} rows)")

    vectors = embed_sequences(
        unique_seqs,
        model_name=model_name,
        batch_size=batch_size,
        window=window,
        device=device,
    )

    seq_to_vec = {seq: vec.tolist() for seq, vec in zip(unique_seqs, vectors)}
    df[EMBEDDING_COLUMN] = df[SEQUENCE_COLUMN].map(seq_to_vec)
    return df


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True, help="Parquet with a `sequence` column")
    ap.add_argument("--out", required=True, help="Parquet with an added `sequence_embedding` column")
    ap.add_argument("--model", default="facebook/esm2_t6_8M_UR50D", help="ESM2 model id on the HF Hub")
    ap.add_argument("--batch-size", type=int, default=32, help="windows per forward pass")
    ap.add_argument("--window", type=int, default=WINDOW_RESIDUES, help="residues per window; long sequences are split and mean-pooled across all windows (no truncation)")
    ap.add_argument("--device", default=None, help="force 'cuda' or 'cpu' (default: auto)")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    in_path = Path(args.inp)
    if not in_path.exists():
        log(f"missing input: {in_path}")
        return 1

    df = run(
        in_path=in_path,
        model_name=args.model,
        batch_size=args.batch_size,
        window=args.window,
        device=args.device,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, compression="zstd", index=False)

    log(f"done: {len(df)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
