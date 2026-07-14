#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import torch
import sys
import logging

from aikisol import Aikisol
from tqdm.auto import tqdm
from pathlib import Path

CHECKPOINT_PATH =  Path(__file__).resolve().parent.parent / "weights" / "aikisol_v2_canonical_147k_full.pt"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

def get_solubilities(sequences: np.ndarray, checkpoint_path=CHECKPOINT_PATH, script_flag=False) -> np.ndarray:
    """
    Predict solubility for an numpy array of protein sequences.
    """
    tqdm.pandas()
    logger = logging.getLogger(__name__)

    if script_flag:
        tqdm_iters = len(sequences)/100
    else:
        tqdm_iters = 1

    model = Aikisol.from_pretrained(checkpoint=checkpoint_path)

    len_seq = len(sequences)
    solubilities = np.empty(len_seq, dtype=np.float32)

    with torch.no_grad():
        for i, sequence in enumerate(tqdm(sequences, desc="Predicting solubility", miniters=tqdm_iters)):
            result = model.predict(sequence)
            solubilities[i] = result.mean_prob[0]

    torch.cuda.empty_cache()

    return solubilities

def main(input_file, output_file, checkpoint_path):
    tqdm.pandas()

    df = pd.read_parquet(input_file, engine="pyarrow")

    df["solubility"] = get_solubilities(df["sequence"].to_numpy(), checkpoint_path=checkpoint_path, script_flag=True)

    df.to_parquet(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--checkpoint", required=True)

    args = parser.parse_args()

    main(args.input, args.output, args.checkpoint)