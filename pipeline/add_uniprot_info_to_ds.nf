#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Attach UniProt `keywords` and `go` to the merged protein dataset.
 *
 * Runs preprocessing/add_uniprot_info_to_ds.py over the merged Parquet: disprot
 * and mobidb rows join to UniProt directly on their accession id, rcsb-pdb rows
 * go through the SIFTS PDB<->UniProt mapping. Row count is preserved.
 *
 * Inputs come from the two download workflows; run them first (both cache, so
 * a repeat run is free):
 *   nextflow run pipeline/download_uniprot.nf
 *   nextflow run pipeline/download_sifts.nf
 *
 * Examples:
 *   nextflow run pipeline/add_uniprot_info_to_ds.nf
 *   nextflow run pipeline/add_uniprot_info_to_ds.nf --merged_filename other.parquet
 */

// Parameter defaults live in nextflow.config (params block).

process addUniprotInfo {
    tag "add_uniprot_info_to_ds"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path merged_dataset
    path uniprot_parquet
    path sifts_mapping

    output:
    path params.uniprot_info_filename

    script:
    """
    python ${params.uniprot_info_script} \\
        --in "${merged_dataset}" \\
        --uniprot "${uniprot_parquet}" \\
        --sifts "${sifts_mapping}" \\
        --out "${params.uniprot_info_filename}"
    """
}

workflow {
    merged_dataset = file("${params.processed_outdir}/${params.merged_filename}")
    uniprot_parquet = file("${params.uniprot_parquet_outdir}/${params.uniprot_parquet_filename}")
    sifts_mapping = file("${params.sifts_outdir}/${params.sifts_filename}")

    addUniprotInfo(merged_dataset, uniprot_parquet, sifts_mapping)
}
