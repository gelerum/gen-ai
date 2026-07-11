#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Build a DisProt sequence dataset:
 *   Uniprot_ID, organism, taxonomy_id, sequence, disorder_mask
 *
 * disorder_mask is a 0/1 string aligned to sequence:
 *   1 = DisProt structural-state disorder
 *   0 = not annotated as disorder
 */

params.disprot_url = params.disprot_url ?: "https://disprot.org/api/v2/download?format=tsv&release=current&term_ontology=IDPO&term_ontology=GO"
params.disprot_outdir = params.disprot_outdir ?: "${launchDir}/data/raw/disprot"
params.disprot_filename = params.disprot_filename ?: "disprot_current_idpo_go.tsv"
params.uniprot_outdir = params.uniprot_outdir ?: "${launchDir}/data/raw/uniprot"
params.dataset_outdir = params.dataset_outdir ?: "${launchDir}/data/processed/disprot"
params.dataset_filename = params.dataset_filename ?: "disprot_sequence_disorder.parquet"
params.uniprot_base_url = params.uniprot_base_url ?: "https://rest.uniprot.org/uniprotkb"
params.disprot_script = params.disprot_script ?: "${launchDir}/preprocessing/build_disprot_dataset.py"
params.limit = params.limit ?: 0
params.max_forks = params.max_forks ?: 6
params.retries = params.retries ?: 3
params.curl_timeout = params.curl_timeout ?: 120

process downloadDisprotTsv {
    tag "disprot"

    input:
    val(download_url)
    val(raw_outdir)
    val(raw_filename)

    output:
    path "disprot.tsv"

    script:
    """
    set -euo pipefail

    mkdir -p "${raw_outdir}"

    curl --fail --location --silent --show-error \\
        --retry ${params.retries} \\
        --max-time ${params.curl_timeout} \\
        --output disprot.tsv \\
        "${download_url}"

    cp disprot.tsv "${raw_outdir}/${raw_filename}"
    """
}

process extractDisprotProteins {
    tag "disprot_proteins"

    input:
    path disprot_tsv

    output:
    path "disprot_proteins.tsv"

    script:
    """
set -euo pipefail

python "${params.disprot_script}" extract-proteins \\
    --disprot-tsv "${disprot_tsv}" \\
    --out disprot_proteins.tsv \\
    --limit ${params.limit}
    """
}

process downloadFastaAndBuildRow {
    tag "${acc}"
    maxForks params.max_forks.toInteger()

    input:
    tuple val(acc), val(organism), val(taxonomy_id), val(starts), val(ends), val(sequence_outdir)

    output:
    path "${acc}.dataset_row.tsv", emit: rows
    path "${acc}.dataset_error.tsv", emit: errors

    script:
    """
set -uo pipefail

mkdir -p "${sequence_outdir}"

fasta="${sequence_outdir}/${acc}.fasta"
fasta_tmp="\${fasta}.tmp.\$\$"
url="${params.uniprot_base_url}/${acc}.fasta"

touch "${acc}.dataset_row.tsv"
touch "${acc}.dataset_error.tsv"

if [ ! -s "\$fasta" ]; then
    curl --fail --location --silent --show-error \\
        --retry ${params.retries} \\
        --max-time ${params.curl_timeout} \\
        --output "\$fasta_tmp" \\
        "\$url"
    curl_status=\$?
    if [ "\$curl_status" -ne 0 ]; then
        rm -f "\$fasta_tmp"
        printf '%s\\t%s\\t%s\\n' "${acc}" "curl_failed" "\$url" > "${acc}.dataset_error.tsv"
        exit 0
    fi
    mv "\$fasta_tmp" "\$fasta"
fi

python "${params.disprot_script}" build-row \\
    --acc "${acc}" \\
    --organism "${organism}" \\
    --taxonomy-id "${taxonomy_id}" \\
    --starts "${starts}" \\
    --ends "${ends}" \\
    --fasta "\$fasta" \\
    --row-out "${acc}.dataset_row.tsv" \\
    --error-out "${acc}.dataset_error.tsv"
    """
}

process writeParquetDataset {
    tag "disprot_parquet"

    input:
    path rows
    path errors
    val(dataset_outdir)
    val(dataset_filename)

    output:
    path "disprot_dataset.manifest.tsv"

    script:
    """
set -euo pipefail

mkdir -p "${dataset_outdir}"

python "${params.disprot_script}" write-parquet \\
    --rows-dir . \\
    --errors-dir . \\
    --out "${dataset_outdir}/${dataset_filename}" \\
    --errors-out "${dataset_outdir}/disprot_sequence_disorder_errors.tsv" \\
    --manifest disprot_dataset.manifest.tsv
    """
}

workflow {
    raw_outdir = params.disprot_outdir.startsWith("/") ? params.disprot_outdir : "${launchDir}/${params.disprot_outdir}"
    sequence_outdir = params.uniprot_outdir.startsWith("/") ? params.uniprot_outdir : "${launchDir}/${params.uniprot_outdir}"
    dataset_outdir = params.dataset_outdir.startsWith("/") ? params.dataset_outdir : "${launchDir}/${params.dataset_outdir}"

    downloadDisprotTsv(params.disprot_url, raw_outdir, params.disprot_filename)

    extractDisprotProteins(downloadDisprotTsv.out)

    protein_rows = extractDisprotProteins.out
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            tuple(
                row["Uniprot_ID"],
                row["organism"],
                row["taxonomy_id"],
                row["starts"],
                row["ends"],
                sequence_outdir
            )
        }

    downloadFastaAndBuildRow(protein_rows)

    downloadFastaAndBuildRow.out.rows
        .collect()
        .set { dataset_rows }

    downloadFastaAndBuildRow.out.errors
        .collect()
        .set { dataset_errors }

    writeParquetDataset(dataset_rows, dataset_errors, dataset_outdir, params.dataset_filename)
}
