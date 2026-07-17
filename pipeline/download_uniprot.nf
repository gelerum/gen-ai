#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download the UniProt/SwissProt release and convert it to Parquet.
 *
 * Fetches the reviewed `uniprot_sprot.dat.gz` flat file into data/raw/uniprot,
 * then runs preprocessing/dat_to_parquet.py over it to publish
 * data/processed/uniprot/uniprot_sprot.parquet with the `keywords` and `go`
 * columns that add_uniprot_info_to_ds.nf joins onto the merged dataset.
 *
 * The download is ~600 MB and the parse takes a while; both are cached by
 * Nextflow, so a re-run with -resume is free.
 *
 * Examples:
 *   nextflow run pipeline/download_uniprot.nf
 *   nextflow run pipeline/download_uniprot.nf --uniprot_dat_url <mirror>
 */

// Parameter defaults live in nextflow.config (params block).

process downloadUniprotDat {
    tag "download_uniprot_dat"
    publishDir params.uniprot_outdir, mode: "copy"

    output:
    path params.uniprot_dat_filename

    script:
    """
    curl -L \\
         --fail \\
         --progress-bar \\
         --retry ${params.retries} \\
         -o "${params.uniprot_dat_filename}" \\
         "${params.uniprot_dat_url}"
    """
}

process uniprotDatToParquet {
    tag "uniprot_dat_to_parquet"
    publishDir params.uniprot_parquet_outdir, mode: "copy"

    input:
    path uniprot_dat

    output:
    path params.uniprot_parquet_filename

    script:
    """
    python ${params.dat_to_parquet_script} \\
        --in "${uniprot_dat}" \\
        --out "${params.uniprot_parquet_filename}"
    """
}

workflow {
    uniprot_dat = downloadUniprotDat()

    uniprotDatToParquet(uniprot_dat)
}
