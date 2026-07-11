#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Stage 2 (RCSB PDB): build the per-protein dataset from the local mmCIF mirror.
 *
 * Runs preprocessing/extract_pdb_features.py over data/raw and publishes a Parquet file
 * (per-protein sequence + disorder mask + CA B-factor) to data/processed.
 * Run this after the PDB download pipeline has populated data/raw.
 *
 * Examples:
 *   nextflow run pipeline/extract_pdb_features.nf
 *   nextflow run pipeline/extract_pdb_features.nf --raw data/raw/pdb_mmCIF --limit 200
 */

// Parameter defaults live in nextflow.config (params block); override with --raw, --limit, etc.

process extractPdbFeatures {
    tag "extract_pdb_features"
    publishDir params.processed_outdir, mode: "copy"

    input:
    val raw_dir

    output:
    path "pdb_protein_features.parquet"

    script:
    def limit_arg = params.limit.toInteger() > 0 ? "--limit ${params.limit}" : ""
    """
    python ${params.extract_script} \\
        --raw "${raw_dir}" \\
        --out pdb_protein_features.parquet \\
        ${limit_arg}
    """
}

workflow {
    raw_dir = params.raw.startsWith('/') ? params.raw : "${launchDir}/${params.raw}"
    extractPdbFeatures(raw_dir)
}
