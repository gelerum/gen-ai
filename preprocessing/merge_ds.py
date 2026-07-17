#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

dataset_path = Path("../data/processed/merged_protein_dataset_human_esm2.parquet")
uniprot_path = Path("../data/processed/uniprot/uniprot_sprot.parquet")
output_path = Path("../data/processed/merged_dataset_with_uniprot.parquet")
sifts_path = "../data/pdb_chain_uniprot.csv.gz"

df = pd.read_parquet(dataset_path)
uniprot = pd.read_parquet(uniprot_path)
sifts = pd.read_csv(sifts_path, sep=",", compression="gzip", comment="#", low_memory=False)

uniprot = uniprot[["accession", "keywords", "go"]].drop_duplicates(subset="accession")

mapping = (sifts[["PDB", "SP_PRIMARY"]].drop_duplicates())
mapping = mapping.rename(columns={"PDB": "id", "SP_PRIMARY": "accession"})

df1 = df[(df['source'] == "mobidb") | (df['source'] == "disprot")]
df2 = df[df['source'] == "rcsb-pdb"]

df2 = df2.merge(mapping, on=["id"], how="left")
merged2 = df2.merge(uniprot, on=["accession"], how="left")
merged2.drop("accession", axis=1, inplace=True)

uniprot = uniprot.rename(columns={"accession": "id"})
merged1 = df1.merge(uniprot, on="id", how="left")

merged = pd.concat([merged1, merged2], ignore_index=True)

output_path.parent.mkdir(parents=True,exist_ok=True,)

merged["transport_terms"] = merged["go"].apply(transport_go_terms)

merged.to_parquet(output_path, compression="zstd", index=False,)