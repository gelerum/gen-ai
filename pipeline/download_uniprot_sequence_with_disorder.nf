#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download a UniProt FASTA sequence and summarize DisProt disorder content.
 *
 * Example:
 *   nextflow run pipeline/download_uniprot_sequence_with_disorder.nf --uniprot_acc P03265
 */

params.uniprot_acc = params.uniprot_acc ?: null
params.disprot_tsv = params.disprot_tsv ?: "${launchDir}/data/raw/disprot/disprot_current_idpo_go.tsv"
params.uniprot_outdir = params.uniprot_outdir ?: "${launchDir}/data/raw/uniprot"
params.disorder_outdir = params.disorder_outdir ?: "${launchDir}/data/processed/disprot"
params.uniprot_base_url = params.uniprot_base_url ?: "https://rest.uniprot.org/uniprotkb"
params.retries = params.retries ?: 3
params.curl_timeout = params.curl_timeout ?: 120

process downloadUniProtSequenceWithDisorder {
    tag "${acc}"

    input:
    val(acc)
    path(disprot_tsv)
    val(sequence_outdir)
    val(summary_outdir)

    output:
    path "${acc}_sequence_disorder.manifest.tsv"

    script:
    """
    set -euo pipefail

    mkdir -p "${sequence_outdir}" "${summary_outdir}"

    fasta="${sequence_outdir}/${acc}.fasta"
    fasta_tmp="\${fasta}.tmp.\$\$"
    summary="${summary_outdir}/${acc}_disorder_summary.tsv"
    summary_tmp="\${summary}.tmp.\$\$"
    url="${params.uniprot_base_url}/${acc}.fasta"

    curl --fail --location --silent --show-error \\
        --retry ${params.retries} \\
        --max-time ${params.curl_timeout} \\
        --output "\$fasta_tmp" \\
        "\$url"

    mv "\$fasta_tmp" "\$fasta"

    fasta_length="\$(awk '
        /^>/ { next }
        { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
        END { print n + 0 }
    ' "\$fasta")"

    awk -F '\\t' -v acc="${acc}" -v fasta_length="\$fasta_length" '
        BEGIN {
            protein_name = "NA"
            gene_name = "NA"
            disprot_id = "NA"
            sequence_length = "NA"
            disorder_content = "NA"
            total_rows = 0
            disorder_regions = 0
        }
        NR == 1 { next }
        \$1 == acc {
            total_rows++
            protein_name = \$3
            gene_name = \$4
            disprot_id = \$2
            sequence_length = \$5
            disorder_content = \$8
            if (\$12 == "Structural state" && \$14 == "disorder") {
                disorder_regions++
            }
        }
        END {
            if (total_rows == 0) {
                printf "No DisProt rows found for UniProt ACC: %s\\n", acc > "/dev/stderr"
                exit 42
            }

            printf "UniProt ACC\\tDisProt ID\\tProtein name\\tGene name\\tSequence length DisProt\\tSequence length FASTA\\tProtein Disorder Content\\tDisorder percent\\tDisProt rows\\tDisorder structural regions\\n"
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%.2f\\t%d\\t%d\\n",
                acc,
                disprot_id,
                protein_name,
                gene_name,
                sequence_length,
                fasta_length,
                disorder_content,
                disorder_content * 100,
                total_rows,
                disorder_regions
        }
    ' "${disprot_tsv}" > "\$summary_tmp"

    mv "\$summary_tmp" "\$summary"

    {
        printf 'source\\tfile\\n'
        printf 'uniprot_fasta\\t%s\\n' "\$fasta"
        printf 'disorder_summary\\t%s\\n' "\$summary"
    } > "${acc}_sequence_disorder.manifest.tsv"
    """
}

workflow {
    if (!params.uniprot_acc) {
        error "Required parameter missing: --uniprot_acc, for example --uniprot_acc P03265"
    }

    acc = params.uniprot_acc.toUpperCase()
    sequence_outdir = params.uniprot_outdir.startsWith('/') ? params.uniprot_outdir : "${launchDir}/${params.uniprot_outdir}"
    summary_outdir = params.disorder_outdir.startsWith('/') ? params.disorder_outdir : "${launchDir}/${params.disorder_outdir}"
    disprot_tsv_path = params.disprot_tsv.startsWith('/') ? params.disprot_tsv : "${launchDir}/${params.disprot_tsv}"
    disprot_tsv = Channel.fromPath(disprot_tsv_path, checkIfExists: true)

    downloadUniProtSequenceWithDisorder(
        acc,
        disprot_tsv,
        sequence_outdir,
        summary_outdir
    )
}
