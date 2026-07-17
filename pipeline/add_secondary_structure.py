#!/usr/bin/env python
"""Add an S4PRED 3-state secondary-structure column to a dataset.

Reads a Parquet file that has a `sequence` column of amino-acid strings, runs
S4PRED over every sequence, and writes the dataset back with an extra column
(default `ss3`) holding a per-residue secondary-structure string aligned 1:1 to
`sequence`:

    H = alpha helix
    E = beta strand
    C = coil / loop
"""

import argparse
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def sanitize(seq) -> str:
    """Upper-case a sequence and replace any non-standard residue with X.

    S4PRED only understands the 20 canonical amino acids; mapping everything
    else to X keeps the output length aligned to the original sequence.
    """
    if not isinstance(seq, str):
        return ""
    seq = seq.strip().upper()
    return "".join(c if c in STANDARD_AA else "X" for c in seq)


def parse_ss_file(path: Path):
    """Extract the secondary-structure string from an S4PRED output file.

    Handles both the `fas` format (header / sequence / ss) and the `ss2`
    vertical format (`  1 M C  0.99 ...`).
    """
    lines = path.read_text().splitlines()

    if any(l.startswith("# PSIPRED") for l in lines):
        ss = []
        for l in lines:
            parts = l.split()
            if len(parts) >= 3 and parts[0].isdigit():
                ss.append(parts[2])
        return "".join(ss)

    nonempty = [l.strip() for l in lines if l.strip()]
    if len(nonempty) >= 3 and nonempty[0].startswith(">"):
        return nonempty[2]
    return None


def predict_batch(seqs, s4pred_dir: Path, device: str):
    """Run S4PRED once over `seqs`; return their SS strings in input order."""
    run_model = s4pred_dir / "run_model.py"

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / "batch.fasta"
        outdir = tmp / "out"
        outdir.mkdir()

        with fasta.open("w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">{i}\n{s}\n")

        cmd = [
            sys.executable, str(run_model),
            "--outfmt", "fas",
            "--save-files", "--save-by-idx", "--silent",
            "--outdir", str(outdir),
            "--device", device,
            str(fasta),
        ]
        subprocess.run(cmd, check=True, cwd=str(s4pred_dir))

        results = []
        for i in range(len(seqs)):
            # --save-by-idx writes s4_out_<i>.<ext>; glob so we don't depend on ext.
            matches = list(outdir.glob(f"s4_out_{i}.*"))
            results.append(parse_ss_file(matches[0]) if matches else None)
        return results


def main(input_file, output_file, s4pred_dir, device, batch_size, column):
    s4pred_dir = Path(s4pred_dir).resolve()

    df = pd.read_parquet(input_file, engine="pyarrow")

    clean = df["sequence"].map(sanitize)
    valid_idx = [i for i, s in zip(df.index, clean) if s]
    seqs = [clean.loc[i] for i in valid_idx]

    logger.info("Predicting secondary structure for %d/%d sequences", len(seqs), len(df))

    if batch_size and batch_size > 0:
        batches = [
            (valid_idx[i:i + batch_size], seqs[i:i + batch_size])
            for i in range(0, len(seqs), batch_size)
        ]
    else:
        batches = [(valid_idx, seqs)]

    ss_by_idx = {}
    done = 0
    for idxs, batch_seqs in batches:
        preds = predict_batch(batch_seqs, s4pred_dir, device)
        ss_by_idx.update(zip(idxs, preds))
        done += len(batch_seqs)
        logger.info("  %d/%d done", done, len(seqs))

    df[column] = df.index.map(lambda i: ss_by_idx.get(i))

    df.to_parquet(output_file, index=False)
    logger.info("Wrote %s with column '%s'", output_file, column)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--s4pred-dir", required=True)
    parser.add_argument("--device", default="cpu", choices=["cpu", "gpu"])
    parser.add_argument("--batch-size", type=int, default=0,
                        help="Sequences per S4PRED run; 0 = all at once.")
    parser.add_argument("--column", default="ss3")

    args = parser.parse_args()

    main(args.input, args.output, args.s4pred_dir, args.device, args.batch_size, args.column)
