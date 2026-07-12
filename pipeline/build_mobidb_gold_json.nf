#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Download MobiDB Gold JSON Lines archive (.mjson.gz) and convert it to JSON.
 *
 * Default output is a compact per-protein JSON dataset with sequence,
 * curated disorder regions, and a sequence-aligned 0/1 disorder mask.
 *
 * Examples:
 *   nextflow run pipeline/build_mobidb_gold_json.nf
 *   nextflow run pipeline/build_mobidb_gold_json.nf --mobidb_mode full
 */

params.mobidb_url = params.mobidb_url ?: "https://protein.bio.unipd.it/shared/mobidb/mobidb_gold_2022_07.mjson.gz"
params.mobidb_outdir = params.mobidb_outdir ?: "${launchDir}/data/raw/mobidb"
params.mobidb_filename = params.mobidb_filename ?: "mobidb_gold_2022_07.mjson.gz"
params.mobidb_dataset_outdir = params.mobidb_dataset_outdir ?: "${launchDir}/data/processed/mobidb"
params.mobidb_dataset_filename = params.mobidb_dataset_filename ?: "mobidb_gold_2022_07.json"
params.mobidb_script = params.mobidb_script ?: "${launchDir}/preprocessing/build_mobidb_gold_json.py"
params.mobidb_mode = params.mobidb_mode ?: "dataset"
params.mobidb_disorder_key = params.mobidb_disorder_key ?: "curated-disorder-merge"
params.mobidb_disorder_variants = params.mobidb_disorder_variants ?: "curated=curated-disorder-merge;homology=homology-disorder-merge;prediction_priority=prediction-disorder-priority;prediction_mobidb_lite=prediction-disorder-mobidb_lite;all_priority=curated-disorder-merge,homology-disorder-merge,prediction-disorder-priority"
params.mobidb_primary_disorder_variant = params.mobidb_primary_disorder_variant ?: "curated"
params.retries = params.retries ?: 3
params.mobidb_curl_timeout = params.mobidb_curl_timeout ?: 0
params.curl_timeout = params.curl_timeout ?: 120

process downloadMobidbGold {
    tag "mobidb_gold"
    debug true

    input:
    val(download_url)
    val(raw_outdir)
    val(raw_filename)
    val(progress_log)

    output:
    path "mobidb_gold.mjson.gz"

    script:
    """
set -euo pipefail

mkdir -p "${raw_outdir}"

target="${raw_outdir}/${raw_filename}"
tmp="\${target}.part"
progress_log="${progress_log}"
: > "\$progress_log"

log_msg() {
    printf '%s %s\\n' "\$(date '+%Y-%m-%d %H:%M:%S')" "\$*" | tee -a "\$progress_log"
}

if [ ! -s "\$target" ]; then
    log_msg "[downloadMobidbGold] downloading MobiDB Gold archive to \$target"
    done_file="\${tmp}.done"
    rm -f "\$done_file"

    (
        last_size="-1"
        while [ ! -f "\$done_file" ]; do
            if [ -f "\$tmp" ]; then
                size=\$(wc -c < "\$tmp" | tr -d ' ')
                if [ "\$size" != "\$last_size" ]; then
                    mib=\$((size / 1024 / 1024))
                    log_msg "[downloadMobidbGold] downloaded \$size bytes (\$mib MiB)"
                    last_size="\$size"
                fi
            else
                log_msg "[downloadMobidbGold] waiting for download to start"
            fi
            sleep 15
        done
    ) &
    monitor_pid=\$!

    cleanup_download_monitor() {
        touch "\$done_file"
        wait "\$monitor_pid" 2>/dev/null || true
        rm -f "\$done_file"
    }
    trap cleanup_download_monitor EXIT

    curl --fail --location --silent --show-error \\
        --retry ${params.retries} \\
        --continue-at - \\
        --max-time ${params.mobidb_curl_timeout} \\
        --output "\$tmp" \\
        "${download_url}"

    cleanup_download_monitor
    trap - EXIT

    mv "\$tmp" "\$target"
    final_size=\$(wc -c < "\$target" | tr -d ' ')
    final_mib=\$((final_size / 1024 / 1024))
    log_msg "[downloadMobidbGold] download complete: \$target (\$final_size bytes, \$final_mib MiB)"
else
    existing_size=\$(wc -c < "\$target" | tr -d ' ')
    existing_mib=\$((existing_size / 1024 / 1024))
    log_msg "[downloadMobidbGold] MobiDB Gold archive already exists: \$target (\$existing_size bytes, \$existing_mib MiB)"
fi

ln -sf "\$target" mobidb_gold.mjson.gz
    """
}

process convertMobidbGoldToJson {
    tag "mobidb_json"
    debug true

    input:
    path mobidb_archive
    val(dataset_outdir)
    val(dataset_filename)
    val(progress_log)

    output:
    path "mobidb_gold_json.manifest.tsv"

    script:
    """
set -euo pipefail

mkdir -p "${dataset_outdir}"

python "${params.mobidb_script}" \\
    --input "${mobidb_archive}" \\
    --out "${dataset_outdir}/${dataset_filename}" \\
    --manifest mobidb_gold_json.manifest.tsv \\
    --mode "${params.mobidb_mode}" \\
    --disorder-key "${params.mobidb_disorder_key}" \\
    --disorder-variants "${params.mobidb_disorder_variants}" \\
    --primary-disorder-variant "${params.mobidb_primary_disorder_variant}" \\
    --progress-log "${progress_log}"
    """
}

workflow {
    raw_outdir = params.mobidb_outdir.startsWith("/") ? params.mobidb_outdir : "${launchDir}/${params.mobidb_outdir}"
    dataset_outdir = params.mobidb_dataset_outdir.startsWith("/") ? params.mobidb_dataset_outdir : "${launchDir}/${params.mobidb_dataset_outdir}"
    progress_log = "${raw_outdir}/mobidb_gold_progress.log"

    downloadMobidbGold(params.mobidb_url, raw_outdir, params.mobidb_filename, progress_log)
    convertMobidbGoldToJson(downloadMobidbGold.out, dataset_outdir, params.mobidb_dataset_filename, progress_log)
}
