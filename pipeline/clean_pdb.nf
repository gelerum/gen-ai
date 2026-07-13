#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Clean the RCSB PDB per-chain features produced by extract_pdb_features.nf.
 *
 * Runs preprocessing/clean_pdb.py over pdb_protein_features.parquet: it keeps
 * well-measured X-ray chains, normalises B-factors by the Wilson B, applies the
 * quality filters, deduplicates sequences per organism and publishes a cleaned
 * Parquet back to data/processed.
 * Run this after extract_pdb_features.nf has finished.
 *
 * Examples:
 *   nextflow run pipeline/clean_pdb.nf
 *   nextflow run pipeline/clean_pdb.nf --pdb_clean_filename clean.parquet
 */

// Parameter defaults live in nextflow.config (params block); override with --pdb_dataset_filename, etc.

process cleanPdb {
    tag "clean_pdb"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path pdb_features

    output:
    path params.pdb_clean_filename

    script:
    """
    python ${params.clean_pdb_script} \\
        --in "${pdb_features}" \\
        --out "${params.pdb_clean_filename}"
    """
}

workflow {
    pdb_features = "${params.processed_outdir}/${params.pdb_dataset_filename}"

    cleanPdb(file(pdb_features))
}
