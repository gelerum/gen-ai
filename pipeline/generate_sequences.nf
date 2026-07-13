#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Balance the clean PDB dataset with random non-folding sequences.
 *
 * Runs preprocessing/generate_sequences.py over pdb_protein_features_clean.parquet:
 * it learns the amino-acid and length distributions from the real chains, samples
 * an equal number of random sequences, drops the B-factor columns and adds a
 * `folds` label (1 = real PDB, 0 = generated), then publishes the result to
 * data/processed.
 * Run this after clean_pdb.nf has finished.
 *
 * The RNG seed is fixed in generate_sequences.py so runs are reproducible;
 * override with --generate_seed to sample a different set.
 *
 * Examples:
 *   nextflow run pipeline/generate_sequences.nf
 *   nextflow run pipeline/generate_sequences.nf --generate_seed 7
 */

// Parameter defaults live in nextflow.config (params block); override with --folds_filename, etc.

process generateSequences {
    tag "generate_sequences"
    publishDir params.processed_outdir, mode: "copy"

    input:
    path pdb_clean

    output:
    path params.folds_filename

    script:
    """
    python ${params.generate_script} \\
        --in "${pdb_clean}" \\
        --out "${params.folds_filename}" \\
        --seed ${params.generate_seed}
    """
}

workflow {
    pdb_clean = "${params.processed_outdir}/${params.pdb_clean_filename}"

    generateSequences(file(pdb_clean))
}
