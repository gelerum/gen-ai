Пайплайн скачивает датасеты из RCSB PDB, DisProt.

По отдельности обрабатывает каждый датасет, вычленяет нужную информацию из каждого, в каждом информация разная.

# Установка Nextflow

```bash
conda create --name nf-env -c bioconda -c conda-forge nextflow pyarrow
source activate nf-env
nextflow info
```

## RCSB PDBx/mmCIF

# Загрузка RCSB PDB в формате PDBx/mmCIF

Проект скачивает структуры RCSB/wwPDB в формате PDBx/mmCIF (`.cif.gz`) и складывает их в `data/raw/pdb_mmCIF`.

Пайплайн находится в `pipeline/download_pdb_mmcif.nf` и использует `curl` (возможны проблемы при загрузке полного датасета размером в 84 гб! TODO: переписать на rsync) внутри задач Nextflow.

Скачать весь доступный архив PDBx/mmCIF:

```bash
nextflow run pipeline/download_pdb_mmcif.nf
```

Для теста можно скачать первые 100 структур:

```bash
nextflow run pipeline/download_pdb_mmcif.nf --limit 100
```

Скачать только структуры из своего списка:

```bash
nextflow run pipeline/download_pdb_mmcif.nf --ids pdb_ids.txt
```

Файл `pdb_ids.txt` может содержать один PDB ID на строку. Также подойдет таблица, где PDB ID находится в первой колонке.

## Параметры

```bash
nextflow run pipeline/download_pdb_mmcif.nf \
  --outdir data/raw/rcsb \
  --chunk_size 500 \
  --max_forks 6
```

Основные параметры:

- `--outdir` — папка для скачанных `.cif.gz`, по умолчанию `data/raw/rcsb`.
- `--ids` — файл со списком PDB ID.
- `--limit` — скачать только первые N структур, удобно для проверки.
- `--chunk_size` — сколько PDB ID обрабатывает одна задача Nextflow.
- `--max_forks` — сколько задач загрузки запускать параллельно.

## Структура данных

Файлы сохраняются по стандартному двухсимвольному shard из архива wwPDB:

```text
data/raw/
  pdb_mmCIF/
    ab/
      1abc.cif.gz
    zz/
      9zzz.cif.gz
```

Например, структура `1abc` будет сохранена как:

```text
data/raw/rcsb/ab/1abc.cif.gz
```

## DisProt Parquet dataset

Собрать итоговый датасет из DisProt:

```bash
nextflow run pipeline/build_disprot_dataset.nf
```

Для быстрой проверки можно ограничить число белков:

```bash
nextflow run pipeline/build_disprot_dataset.nf --limit 100
```

Пайплайн делает все шаги сам:

1. Скачивает текущий DisProt TSV.
2. Берет уникальные `UniProt ACC`.
3. Скачивает FASTA-последовательности из UniProt.
4. Строит бинарную маску disorder.
5. Сохраняет итоговый датасет в Parquet.

Сырые данные сохраняются сюда:

```text
data/raw/disprot/disprot_current_idpo_go.tsv
data/raw/uniprot/
```

По умолчанию используется DisProt URL:

```text
https://disprot.org/api/v2/download?format=tsv&release=current&term_ontology=IDPO&term_ontology=GO
```

Итоговый файл:

```text
data/processed/disprot/disprot_sequence_disorder.parquet
```

Если для части UniProt ID не удалось получить последовательность или FASTA оказался пустым, эти белки пропускаются, а причины сохраняются отдельно:

```text
data/processed/disprot/disprot_sequence_disorder_errors.tsv
```

Колонки итогового Parquet:

```text
Uniprot_ID
organism
taxonomy_id
sequence
disorder_mask
```

`disorder_mask` — строка той же длины, что и `sequence`: `1` означает disorder-позицию, `0` означает order/не размечено как disorder.
