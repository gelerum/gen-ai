#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.weights_dir = "${launchDir}/weights"

process download_weights {

    publishDir params.weights_dir, mode: 'copy'

    output:
    path "aikisol_v2_canonical_147k_full.pt"

    script:
    def model_name = "aikisol_v2_canonical_147k_full.pt"
    def url = "https://zenodo.org/records/20171586/files/apache_tier__aikisol_v2_canonical_147k_full.pt?download=1"

    """
    curl -L \
         --fail \
         --progress-bar \
         -o ${model_name} \
         "${url}"
    """
}

workflow {
    download_weights()
}