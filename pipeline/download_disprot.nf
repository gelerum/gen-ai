#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download current DisProt annotations as TSV.
 *
 * Example:
 *   nextflow run pipeline/download_disprot.nf
 */

params.disprot_outdir = params.disprot_outdir ?: "${launchDir}/data/raw/disprot"
params.disprot_url = params.disprot_url ?: "https://disprot.org/api/v2/download?format=tsv&release=current&term_ontology=IDPO&term_ontology=GO"
params.disprot_filename = params.disprot_filename ?: "disprot_current_idpo_go.tsv"
params.retries = params.retries ?: 3
params.curl_timeout = params.curl_timeout ?: 120

process downloadDisprotTsv {
    tag "disprot"

    input:
    val(download_url)
    val(download_outdir)
    val(output_name)

    output:
    path "disprot_download.manifest.tsv"

    script:
    """
    set -euo pipefail

    mkdir -p "${download_outdir}"

    target="${download_outdir}/${output_name}"
    tmp="\${target}.tmp.\$\$"

    curl --fail --location --silent --show-error \\
        --retry ${params.retries} \\
        --max-time ${params.curl_timeout} \\
        --output "\$tmp" \\
        "${download_url}"

    mv "\$tmp" "\$target"

    if command -v sha256sum >/dev/null 2>&1; then
        checksum="\$(sha256sum "\$target" | awk '{ print \$1 }')"
    else
        checksum="\$(shasum -a 256 "\$target" | awk '{ print \$1 }')"
    fi

    {
        printf 'source\\turl\\tfile\\tsha256\\n'
        printf 'disprot\\t%s\\t%s\\t%s\\n' "${download_url}" "\$target" "\$checksum"
    } > disprot_download.manifest.tsv
    """
}

workflow {
    download_outdir = params.disprot_outdir.startsWith('/') ? params.disprot_outdir : "${launchDir}/${params.disprot_outdir}"

    downloadDisprotTsv(
        params.disprot_url,
        download_outdir,
        params.disprot_filename
    )
}
