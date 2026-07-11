#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Build a dataset with full UniProt sequences and DisProt disorder percent.
 *
 * Example:
 *   nextflow run pipeline/build_disprot_sequence_disorder_dataset.nf
 *   nextflow run pipeline/build_disprot_sequence_disorder_dataset.nf --limit 100
 */

params.disprot_tsv = params.disprot_tsv ?: "${launchDir}/data/raw/disprot/disprot_current_idpo_go.tsv"
params.uniprot_outdir = params.uniprot_outdir ?: "${launchDir}/data/raw/uniprot"
params.dataset_outdir = params.dataset_outdir ?: "${launchDir}/data/processed/disprot"
params.dataset_filename = params.dataset_filename ?: "sequence_disorder_dataset.tsv"
params.uniprot_base_url = params.uniprot_base_url ?: "https://rest.uniprot.org/uniprotkb"
params.limit = params.limit ?: 0
params.max_forks = params.max_forks ?: 6
params.retries = params.retries ?: 3
params.curl_timeout = params.curl_timeout ?: 120

process extractDisprotProteins {
    tag "disprot_proteins"

    input:
    path disprot_tsv

    output:
    path "disprot_proteins.tsv"

    script:
    """
    set -euo pipefail

    awk -F '\\t' '
        NR == 1 { next }
        {
            acc = \$1
            if (acc == "") next

            rows[acc]++

            if (!(acc in seen)) {
                seen[acc] = 1
                disprot_id[acc] = \$2
                protein_name[acc] = \$3
                gene_name[acc] = \$4
                sequence_length[acc] = \$5
                disorder_content[acc] = \$8
            }

            if (\$12 == "Structural state" && \$14 == "disorder") {
                disorder_regions[acc]++
            }
        }
        END {
            for (acc in seen) {
                disorder_percent = disorder_content[acc] * 100
                printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%.2f\\t%d\\t%d\\n",
                    acc,
                    disprot_id[acc],
                    protein_name[acc],
                    gene_name[acc],
                    sequence_length[acc],
                    disorder_content[acc],
                    disorder_percent,
                    rows[acc],
                    disorder_regions[acc] + 0
            }
        }
    ' "$disprot_tsv" | sort -k1,1 > proteins.body.tsv

    {
        printf 'UniProt ACC\\tDisProt ID\\tProtein name\\tGene name\\tSequence length DisProt\\tProtein Disorder Content\\tDisorder percent\\tDisProt rows\\tDisorder structural regions\\n'
        if [ ${params.limit} -gt 0 ]; then
            head -n ${params.limit} proteins.body.tsv
        else
            cat proteins.body.tsv
        fi
    } > disprot_proteins.tsv
    """
}

process downloadSequenceDatasetRow {
    tag "${acc}"
    maxForks params.max_forks.toInteger()

    input:
    tuple val(acc), val(disprot_id), val(protein_name), val(gene_name),
          val(disprot_length), val(disorder_content), val(disorder_percent),
          val(disprot_rows), val(disorder_regions), val(sequence_outdir)

    output:
    path "${acc}.sequence_disorder.row.tsv"

    script:
    """
    set -euo pipefail

    mkdir -p "${sequence_outdir}"

    fasta="${sequence_outdir}/${acc}.fasta"
    fasta_tmp="\${fasta}.tmp.\$\$"
    url="${params.uniprot_base_url}/${acc}.fasta"

    if [ ! -s "\$fasta" ]; then
        curl --fail --location --silent --show-error \\
            --retry ${params.retries} \\
            --max-time ${params.curl_timeout} \\
            --output "\$fasta_tmp" \\
            "\$url"

        mv "\$fasta_tmp" "\$fasta"
    fi

    sequence="\$(awk '
        /^>/ { next }
        { gsub(/[[:space:]]/, "", \$0); printf "%s", \$0 }
    ' "\$fasta")"

    fasta_length="\${#sequence}"

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${acc}" \\
        "${disprot_id}" \\
        "${protein_name}" \\
        "${gene_name}" \\
        "${disprot_length}" \\
        "\$fasta_length" \\
        "${disorder_content}" \\
        "${disorder_percent}" \\
        "${disprot_rows}" \\
        "${disorder_regions}" \\
        "\$sequence" \\
        > "${acc}.sequence_disorder.row.tsv"
    """
}

process mergeSequenceDisorderDataset {
    tag "sequence_disorder_dataset"

    input:
    path rows
    val(dataset_outdir)
    val(dataset_filename)

    output:
    path "sequence_disorder_dataset.manifest.tsv"

    script:
    """
    set -euo pipefail

    mkdir -p "${dataset_outdir}"

    dataset="${dataset_outdir}/${dataset_filename}"
    dataset_tmp="\${dataset}.tmp.\$\$"

    printf 'UniProt ACC\\tDisProt ID\\tProtein name\\tGene name\\tSequence length DisProt\\tSequence length FASTA\\tProtein Disorder Content\\tDisorder percent\\tDisProt rows\\tDisorder structural regions\\tSequence\\n' > "\$dataset_tmp"
    cat ${rows} | sort -k1,1 >> "\$dataset_tmp"

    mv "\$dataset_tmp" "\$dataset"

    {
        printf 'source\\tfile\\n'
        printf 'sequence_disorder_dataset\\t%s\\n' "\$dataset"
    } > sequence_disorder_dataset.manifest.tsv
    """
}

workflow {
    disprot_tsv_path = params.disprot_tsv.startsWith('/') ? params.disprot_tsv : "${launchDir}/${params.disprot_tsv}"
    sequence_outdir = params.uniprot_outdir.startsWith('/') ? params.uniprot_outdir : "${launchDir}/${params.uniprot_outdir}"
    dataset_outdir = params.dataset_outdir.startsWith('/') ? params.dataset_outdir : "${launchDir}/${params.dataset_outdir}"
    disprot_tsv = Channel.fromPath(disprot_tsv_path, checkIfExists: true)

    extractDisprotProteins(disprot_tsv)

    protein_rows = extractDisprotProteins.out
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row['UniProt ACC'],
                row['DisProt ID'],
                row['Protein name'],
                row['Gene name'],
                row['Sequence length DisProt'],
                row['Protein Disorder Content'],
                row['Disorder percent'],
                row['DisProt rows'],
                row['Disorder structural regions'],
                sequence_outdir
            )
        }

    downloadSequenceDatasetRow(protein_rows)

    downloadSequenceDatasetRow.out
        .collect()
        .set { dataset_rows }

    mergeSequenceDisorderDataset(
        dataset_rows,
        dataset_outdir,
        params.dataset_filename
    )
}
