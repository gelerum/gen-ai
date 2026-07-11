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

python - <<'PY'
import csv

proteins = {}

with open("${disprot_tsv}", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\\t")
    for row in reader:
        acc = row["UniProt ACC"].strip()
        if not acc:
            continue

        item = proteins.setdefault(acc, {
            "Uniprot_ID": acc,
            "organism": row["Organism"],
            "taxonomy_id": row["NCBI Taxon ID"],
            "starts": [],
            "ends": [],
        })

        is_disorder = (
            row["Term namespace"] == "Structural state"
            and row["Term name"] == "disorder"
        )
        if is_disorder:
            try:
                start = int(row["Start"])
                end = int(row["End"])
            except ValueError:
                continue
            if start > 0 and end >= start:
                item["starts"].append(start)
                item["ends"].append(end)

rows = sorted(proteins.values(), key=lambda x: x["Uniprot_ID"])
limit = int("${params.limit}")
if limit > 0:
    rows = rows[:limit]

with open("disprot_proteins.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["Uniprot_ID", "organism", "taxonomy_id", "starts", "ends"],
        delimiter="\\t",
        lineterminator="\\n",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({
            "Uniprot_ID": row["Uniprot_ID"],
            "organism": row["organism"],
            "taxonomy_id": row["taxonomy_id"],
            "starts": ",".join(map(str, row["starts"])),
            "ends": ",".join(map(str, row["ends"])),
        })
PY
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

python - <<'PY'
import csv

acc = "${acc}"
organism = "${organism}"
taxonomy_id = "${taxonomy_id}"
starts = [int(x) for x in "${starts}".split(",") if x and x != "null"]
ends = [int(x) for x in "${ends}".split(",") if x and x != "null"]

sequence_parts = []
with open("${sequence_outdir}/${acc}.fasta") as handle:
    for line in handle:
        line = line.strip()
        if line and not line.startswith(">"):
            sequence_parts.append(line)

sequence = "".join(sequence_parts)
if not sequence:
    with open(f"{acc}.dataset_error.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
        writer.writerow([acc, "empty_fasta", f"${sequence_outdir}/{acc}.fasta"])
    raise SystemExit(0)

mask = ["0"] * len(sequence)

for start, end in zip(starts, ends):
    left = max(start, 1) - 1
    right = min(end, len(sequence))
    for idx in range(left, right):
        mask[idx] = "1"

with open(f"{acc}.dataset_row.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow([acc, organism, taxonomy_id, sequence, "".join(mask)])
PY
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

python - <<'PY'
import csv
import sys
from pathlib import Path

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:
    sys.stderr.write(
        "pyarrow is required to write Parquet. Install it in the Nextflow environment, "
        "for example: conda install -c conda-forge pyarrow\\n"
    )
    raise

rows = []
row_files = sorted(Path(".").glob("*.dataset_row.tsv"))
for row_file in row_files:
    with open(row_file, newline="") as handle:
        reader = csv.reader(handle, delimiter="\\t")
        for values in reader:
            if not values:
                continue
            rows.append({
                "Uniprot_ID": values[0],
                "organism": values[1],
                "taxonomy_id": values[2],
                "sequence": values[3],
                "disorder_mask": values[4],
            })

rows.sort(key=lambda item: item["Uniprot_ID"])

table = pa.Table.from_pylist(
    rows,
    schema=pa.schema([
        ("Uniprot_ID", pa.string()),
        ("organism", pa.string()),
        ("taxonomy_id", pa.string()),
        ("sequence", pa.string()),
        ("disorder_mask", pa.string()),
    ]),
)

out_path = Path("${dataset_outdir}") / "${dataset_filename}"
pq.write_table(table, out_path)

errors = []
error_files = sorted(Path(".").glob("*.dataset_error.tsv"))
for error_file in error_files:
    with open(error_file, newline="") as handle:
        reader = csv.reader(handle, delimiter="\\t")
        for values in reader:
            if not values:
                continue
            errors.append({
                "Uniprot_ID": values[0],
                "error_type": values[1],
                "detail": values[2] if len(values) > 2 else "",
            })

error_path = Path("${dataset_outdir}") / "disprot_sequence_disorder_errors.tsv"
with open(error_path, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["Uniprot_ID", "error_type", "detail"],
        delimiter="\\t",
        lineterminator="\\n",
    )
    writer.writeheader()
    writer.writerows(errors)

with open("disprot_dataset.manifest.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["source", "file", "rows"])
    writer.writerow(["disprot_parquet", str(out_path), len(rows)])
    writer.writerow(["disprot_errors", str(error_path), len(errors)])
PY
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
