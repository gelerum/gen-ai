#!/usr/bin/env python3
"""Balance the clean PDB dataset with random ("non-folding") sequences.

The cleaned PDB dataset holds real, crystallised protein chains — sequences that
demonstrably fold. To train a foldability classifier we need negatives too, so
this step learns two distributions from the real data and samples matching random
sequences from them:

  1. amino-acid weights  — the empirical frequency of each residue letter across
     all real sequences (the "composition" of real proteins)
  2. length distribution — the best fit (lowest AIC) among norm/laplace/expon/
     gamma/lognorm to the real chain lengths

It then generates as many random sequences as there are real rows: draw a length
from (2), fill it with residues sampled i.i.d. from (1). The two sets are stacked
and labelled with a `folds` column (1 = real PDB, 0 = generated). B-factor
columns are dropped here — the negatives have no structure, so they are not part
of this dataset.

Requires: pandas, numpy, scipy, pyarrow.

Examples:
    python preprocessing/generate_sequences.py \\
        --in data/processed/pdb_protein_features_clean.parquet \\
        --out data/processed/pdb_folds_dataset.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as st

# Fixed RNG seed so the generated sequences are reproducible run to run.
DEFAULT_SEED = 42

# Columns holding B-factor data; dropped because generated rows have no structure.
BFACTOR_COLUMNS = ["bfactor", "bfactor_norm"]

# Candidate length distributions, fitted and ranked by AIC (see fit_length_dist).
LENGTH_CANDIDATES = [st.norm, st.laplace, st.expon, st.gamma, st.lognorm]


def log(message: str) -> None:
    print(message, file=sys.stderr)


def amino_acid_weights(sequences: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    """Empirical residue alphabet and its probability weights over all sequences."""
    counts = pd.Series(list("".join(sequences))).value_counts()
    alphabet = counts.index.to_numpy()
    weights = counts.to_numpy(dtype=float)
    weights /= weights.sum()
    log(f"amino-acid alphabet: {len(alphabet)} letters over {counts.sum()} residues")
    return alphabet, weights


def fit_length_dist(lengths: pd.Series):
    """Fit each candidate distribution and return the best (dist, params) by AIC."""
    best = None
    for dist in LENGTH_CANDIDATES:
        params = dist.fit(lengths)
        loglik = np.sum(dist.logpdf(lengths, *params))
        aic = 2 * len(params) - 2 * loglik
        log(f"length fit {dist.name}: AIC={aic:.1f}")
        if best is None or aic < best[0]:
            best = (aic, dist, params)
    _, dist, params = best
    log(f"best length distribution: {dist.name}")
    return dist, params


def generate_sequences(
    n: int,
    alphabet: np.ndarray,
    weights: np.ndarray,
    length_dist,
    length_params,
    rng: np.random.Generator,
) -> list[str]:
    """Sample `n` random sequences: lengths from the fit, residues from the weights."""
    lengths = length_dist.rvs(*length_params, size=n, random_state=rng)
    lengths = np.clip(np.round(lengths).astype(int), 1, None)

    # Draw every residue in one call, then slice it into per-sequence strings.
    residues = rng.choice(alphabet, size=int(lengths.sum()), p=weights)
    sequences = []
    start = 0
    for length in lengths:
        sequences.append("".join(residues[start : start + length]))
        start += length
    return sequences


def build(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    log(f"loaded: {len(df)} rows")

    # Real, folding proteins: drop B-factor data and label them positive.
    real = df.drop(columns=[c for c in BFACTOR_COLUMNS if c in df.columns])
    real = real.copy()
    real["folds"] = 1

    # Learn the distributions from the real data, then sample matching negatives.
    alphabet, weights = amino_acid_weights(real["sequence"].dropna().astype(str))
    length_dist, length_params = fit_length_dist(real["length"].dropna())

    sequences = generate_sequences(
        n=len(real),
        alphabet=alphabet,
        weights=weights,
        length_dist=length_dist,
        length_params=length_params,
        rng=rng,
    )
    generated = pd.DataFrame(
        {
            "sequence": sequences,
            "length": [len(s) for s in sequences],
            "folds": 0,
        }
    )
    log(f"generated {len(generated)} random sequences")

    # Stack; metadata columns absent from `generated` become NaN for those rows.
    combined = pd.concat([real, generated], ignore_index=True)
    log(f"combined: {len(combined)} rows ({int(combined['folds'].sum())} folds=1, "
        f"{int((combined['folds'] == 0).sum())} folds=0)")
    return combined


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True, help="clean PDB Parquet to read")
    ap.add_argument("--out", required=True, help="Parquet with a `folds` label column")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED, help=f"RNG seed (default: {DEFAULT_SEED})")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    in_path = Path(args.inp)
    if not in_path.exists():
        log(f"missing input: {in_path}")
        return 1

    rng = np.random.default_rng(args.seed)
    df = build(pd.read_parquet(in_path), rng)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, compression="zstd", index=False)

    log(f"done: {len(df)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
