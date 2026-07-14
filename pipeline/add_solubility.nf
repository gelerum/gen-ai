nextflow.enable.dsl = 2

params.input      = "${launchDir}/data/processed/merged_protein_dataset_clean.parquet"
params.checkpoint = "${launchDir}/weights/aikisol_v2_canonical_147k_full.pt"

process compute_solubility {

    publishDir "${launchDir}/data/processed", mode: 'copy'

    input:
    path input_parquet
    path checkpoint

    output:
    path "ds_solubility.parquet"

    script:
    """
    python ${launchDir}/pipeline/add_solubility.py \
        --input $input_parquet \
        --output ds_solubility.parquet \
        --checkpoint $checkpoint
    """
}

workflow {
    input_ch = Channel.fromPath(params.input)
    checkpoint_ch = Channel.fromPath(params.checkpoint)

    compute_solubility(input_ch, checkpoint_ch)
}