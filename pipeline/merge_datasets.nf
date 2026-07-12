#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Merge the three processed datasets (DisProt, MobiDB, RCSB PDB) into one
 * unified per-protein Parquet file with a `source` column.
 *
 * Runs preprocessing/merge_datasets.py over the already-built Parquet datasets
 * in data/processed and publishes the merged file back to data/processed.
 * Run this after the DisProt, MobiDB, and PDB feature pipelines have finished.
 *
 * Examples:
 *   nextflow run pipeline/merge_datasets.nf
 *   nextflow run pipeline/merge_datasets.nf --merged_filename merged.parquet
 */

// Parameter defaults live in nextflow.config (params block); override with --disprot_dataset, etc.

process mergeDatasets {
    tag "merge_datasets"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path disprot_dataset
    path mobidb_dataset
    path pdb_dataset

    output:
    path params.merged_filename

    script:
    """
    python ${params.merge_script} \\
        --disprot "${disprot_dataset}" \\
        --mobidb "${mobidb_dataset}" \\
        --pdb "${pdb_dataset}" \\
        --out "${params.merged_filename}"
    """
}

workflow {
    disprot_dataset = "${params.dataset_outdir}/${params.dataset_filename}"
    mobidb_dataset = "${params.mobidb_dataset_outdir}/${params.mobidb_dataset_filename}"
    pdb_dataset = "${params.processed_outdir}/${params.pdb_dataset_filename}"

    mergeDatasets(
        file(disprot_dataset),
        file(mobidb_dataset),
        file(pdb_dataset),
    )
}
