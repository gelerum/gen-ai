#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Attach per-sequence ESM2 embeddings to one organism of the clean dataset.
 *
 * Runs preprocessing/embed_esm2.py over the cleaned merged Parquet: it keeps
 * the rows for a single taxonomy_id (9606 = human by default) and adds one
 * mean-pooled ESM2 embedding vector per protein sequence, then publishes the
 * result back to data/processed.
 * Run this after clean_merged_dataset.nf has finished.
 *
 * The ESM2 model is pulled from the Hugging Face Hub the first time this stage
 * runs (nothing is downloaded until then); it uses CUDA when available and
 * falls back to CPU otherwise.
 *
 * Examples:
 *   nextflow run pipeline/embed_esm2.nf
 *   nextflow run pipeline/embed_esm2.nf --embed_taxonomy_id 10090
 */

// Parameter defaults live in nextflow.config (params block); override with --embed_model, etc.

process embedEsm2 {
    tag "embed_esm2"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path clean_dataset

    output:
    path params.embed_filename

    script:
    """
    python ${params.embed_script} \\
        --in "${clean_dataset}" \\
        --out "${params.embed_filename}" \\
        --taxonomy-id "${params.embed_taxonomy_id}" \\
        --model "${params.embed_model}" \\
        --batch-size ${params.embed_batch_size} \\
        --window ${params.embed_window}
    """
}

workflow {
    clean_dataset = "${params.processed_outdir}/${params.clean_filename}"

    embedEsm2(file(clean_dataset))
}
