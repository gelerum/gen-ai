#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Clean the unified merged dataset produced by merge_datasets.nf.
 *
 * Runs preprocessing/clean_merged_dataset.py over the merged Parquet: it drops
 * empty/label-less/duplicate rows and canonicalises organism/sequence casing,
 * then publishes a cleaned Parquet back to data/processed.
 * Run this after merge_datasets.nf has finished.
 *
 * Examples:
 *   nextflow run pipeline/clean_merged_dataset.nf
 *   nextflow run pipeline/clean_merged_dataset.nf --clean_filename clean.parquet
 */

// Parameter defaults live in nextflow.config (params block); override with --merged_filename, etc.

process cleanMergedDataset {
    tag "clean_merged_dataset"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path merged_dataset

    output:
    path params.clean_filename

    script:
    """
    python ${params.clean_script} \\
        --in "${merged_dataset}" \\
        --out "${params.clean_filename}"
    """
}

workflow {
    merged_dataset = "${params.processed_outdir}/${params.merged_filename}"

    cleanMergedDataset(file(merged_dataset))
}
