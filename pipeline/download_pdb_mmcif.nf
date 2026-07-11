#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download RCSB/wwPDB structures in PDBx/mmCIF format (.cif.gz).
 *
 * Examples:
 *   nextflow run pipeline/download_pdb_mmcif.nf
 *   nextflow run pipeline/download_pdb_mmcif.nf --limit 100
 *   nextflow run pipeline/download_pdb_mmcif.nf --ids pdb_ids.txt --outdir data/raw
 */

params.outdir = params.outdir ?: "${launchDir}/data/raw"
params.ids = params.ids ?: null
params.limit = params.limit ?: 0
params.chunk_size = params.chunk_size ?: 200
params.max_forks = params.max_forks ?: 6
params.retries = params.retries ?: 3
params.curl_timeout = params.curl_timeout ?: 120
params.pdb_index_url = params.pdb_index_url ?: "https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
params.mmcif_base_url = params.mmcif_base_url ?: "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF"

process fetchPdbIndex {
    tag "entries.idx"

    output:
    path "entries.idx"

    script:
    """
    curl --fail --location --silent --show-error \\
        --retry ${params.retries} \\
        --max-time ${params.curl_timeout} \\
        --output entries.idx \\
        "${params.pdb_index_url}"
    """
}

process normalizePdbIds {
    tag "pdb_ids"

    input:
    path source_ids

    output:
    path "pdb_ids.txt"

    script:
    """
    awk '
        /^[[:space:]]*\$/ { next }
        /^#/ { next }
        /^IDCODE/ { next }
        /^-/ { next }
        {
            id = tolower(\$1)
            if (id ~ /^[0-9][a-z0-9]{3}\$/) print id
        }
    ' "$source_ids" | sort -u > all_ids.txt

    if [ ${params.limit} -gt 0 ]; then
        head -n ${params.limit} all_ids.txt > pdb_ids.txt
    else
        mv all_ids.txt pdb_ids.txt
    fi
    """
}

process downloadPdbMmcif {
    tag "chunk_${chunk_id}"
    maxForks params.max_forks.toInteger()

    input:
    tuple val(chunk_id), path(ids_file), val(download_outdir)

    output:
    path "chunk_${chunk_id}.manifest.tsv"

    script:
    """
    set -euo pipefail

    outdir="${download_outdir}"
    base_url="${params.mmcif_base_url}"
    mkdir -p "\$outdir"
    : > "chunk_${chunk_id}.manifest.tsv"

    while IFS= read -r pdb_id; do
        [ -n "\$pdb_id" ] || continue

        shard="\${pdb_id:1:2}"
        target_dir="\$outdir/\$shard"
        target="\$target_dir/\${pdb_id}.cif.gz"
        tmp="\${target}.tmp.\$\$"
        url="\$base_url/\$shard/\${pdb_id}.cif.gz"

        mkdir -p "\$target_dir"

        if [ -s "\$target" ]; then
            printf '%s\\t%s\\t%s\\n' "\$pdb_id" "exists" "\$target" >> "chunk_${chunk_id}.manifest.tsv"
            continue
        fi

        curl --fail --location --silent --show-error \\
            --retry ${params.retries} \\
            --continue-at - \\
            --max-time ${params.curl_timeout} \\
            --output "\$tmp" \\
            "\$url"

        mv "\$tmp" "\$target"
        printf '%s\\t%s\\t%s\\n' "\$pdb_id" "downloaded" "\$target" >> "chunk_${chunk_id}.manifest.tsv"
    done < "$ids_file"
    """
}

workflow {
    download_outdir = params.outdir.startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"

    if (params.ids) {
        source_ids = Channel.fromPath(params.ids, checkIfExists: true)
    } else {
        fetchPdbIndex()
        source_ids = fetchPdbIndex.out
    }

    normalizePdbIds(source_ids)

    normalizePdbIds.out
        .splitText(by: params.chunk_size as int, file: true)
        .ifEmpty { error "No PDB IDs found" }
        .map { chunk_file -> tuple(chunk_file.baseName.replaceAll(/\D+/, '') ?: chunk_file.name, chunk_file, download_outdir) }
        | downloadPdbMmcif
}
