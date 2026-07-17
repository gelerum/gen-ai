#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download the S4PRED secondary-structure predictor.
 *
 * Clones the S4PRED repo and fetches its model weights (~430MB uncompressed),
 * then publishes the ready-to-run directory to tools/s4pred.
 * Run add_secondary_structure.nf afterwards to use it.
 *
 * S4PRED itself needs PyTorch and BioPython at prediction time.
 *
 * Example:
 *   nextflow run pipeline/download_s4pred.nf
 */

params.tools_dir          = "${launchDir}/tools"
params.s4pred_repo        = "https://github.com/psipred/s4pred"
params.s4pred_weights_url = "http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights.tar.gz"

process download_s4pred {

    publishDir params.tools_dir, mode: 'copy'

    output:
    path "s4pred"

    script:
    """
    set -euo pipefail

    git clone --depth 1 ${params.s4pred_repo} s4pred

    curl -L --fail --progress-bar \
         -o s4pred/weights.tar.gz \
         "${params.s4pred_weights_url}"

    tar -xzf s4pred/weights.tar.gz -C s4pred
    rm -f s4pred/weights.tar.gz
    """
}

workflow {
    download_s4pred()
}
