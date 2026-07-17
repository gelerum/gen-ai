#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download the SIFTS PDB <-> UniProt residue-level mapping.
 *
 * Fetches `pdb_chain_uniprot.csv.gz` from the EBI SIFTS flatfiles into
 * data/raw/sifts. add_uniprot_info_to_ds.nf uses it to resolve the PDB entry
 * ids of the rcsb-pdb rows to UniProt accessions; the disprot and mobidb rows
 * already carry accessions as their id and need no mapping.
 *
 * This file used to be committed into git as a 6 MB blob; it is raw input data
 * and belongs in data/raw, which is gitignored.
 *
 * Examples:
 *   nextflow run pipeline/download_sifts.nf
 */

// Parameter defaults live in nextflow.config (params block).

process downloadSifts {
    tag "download_sifts"
    publishDir params.sifts_outdir, mode: "copy"

    output:
    path params.sifts_filename

    script:
    """
    curl -L \\
         --fail \\
         --progress-bar \\
         --retry ${params.retries} \\
         -o "${params.sifts_filename}" \\
         "${params.sifts_url}"
    """
}

workflow {
    downloadSifts()
}
