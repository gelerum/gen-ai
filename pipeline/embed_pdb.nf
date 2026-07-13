#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Attach per-sequence ESM2 embeddings to the cleaned RCSB PDB dataset.
 *
 * Runs preprocessing/embed_sequences.py over pdb_protein_features_clean.parquet:
 * it adds one mean-pooled ESM2 `sequence_embedding` vector per row and publishes
 * the result back to data/processed.
 * Run this after clean_pdb.nf has finished.
 *
 * The ESM2 model is pulled from the Hugging Face Hub the first time this stage
 * runs (nothing is downloaded until then); it uses CUDA when available and
 * falls back to CPU otherwise.
 *
 * Examples:
 *   nextflow run pipeline/embed_pdb.nf
 *   nextflow run pipeline/embed_pdb.nf --embed_model facebook/esm2_t12_35M_UR50D
 */

// Parameter defaults live in nextflow.config (params block); override with --embed_model, etc.

process embedPdb {
    tag "embed_pdb"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path pdb_clean

    output:
    path params.embed_pdb_filename

    script:
    """
    python ${params.embed_sequences_script} \\
        --in "${pdb_clean}" \\
        --out "${params.embed_pdb_filename}" \\
        --model "${params.embed_model}" \\
        --batch-size ${params.embed_batch_size} \\
        --window ${params.embed_window}
    """
}

workflow {
    pdb_clean = "${params.processed_outdir}/${params.pdb_clean_filename}"

    embedPdb(file(pdb_clean))
}
