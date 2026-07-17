#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Add an S4PRED 3-state secondary-structure column to a dataset.
 *
 * Takes any Parquet with a `sequence` column and appends `ss3`: a per-residue
 * H/E/C string (helix / strand / coil) aligned to each sequence.
 * Run download_s4pred.nf first so tools/s4pred exists.
 *
 * Examples:
 *   nextflow run pipeline/add_secondary_structure.nf
 *   nextflow run pipeline/add_secondary_structure.nf --ss_device gpu
 */

params.input         = "${launchDir}/data/processed/merged_protein_dataset_clean.parquet"
params.s4pred_dir    = "${launchDir}/tools/s4pred"
params.ss_device     = "cpu"
params.ss_batch_size = 0
params.ss_column     = "ss3"

process add_secondary_structure {

    publishDir "${launchDir}/data/processed", mode: 'copy'

    input:
    path input_parquet
    path s4pred_dir

    output:
    path "ds_secondary_structure.parquet"

    script:
    """
    python ${launchDir}/pipeline/add_secondary_structure.py \
        --input $input_parquet \
        --output ds_secondary_structure.parquet \
        --s4pred-dir $s4pred_dir \
        --device ${params.ss_device} \
        --batch-size ${params.ss_batch_size} \
        --column ${params.ss_column}
    """
}

workflow {
    input_ch  = Channel.fromPath(params.input)
    s4pred_ch = Channel.fromPath(params.s4pred_dir, type: 'dir')

    add_secondary_structure(input_ch, s4pred_ch)
}
