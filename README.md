# Загрузка сырых биологических данных

Проект скачивает сырые данные из внешних источников и складывает их в `data/raw`.

Сейчас есть несколько пайплайнов:

- `pipeline/download_pdb_mmcif.nf` — структуры RCSB/wwPDB в формате PDBx/mmCIF (`.cif.gz`).
- `pipeline/download_disprot.nf` — текущий TSV export DisProt.
- `pipeline/build_disprot_sequence_disorder_dataset.nf` — датасет `[последовательность, disorder процент]` для белков из DisProt.
- `pipeline/download_uniprot_sequence_with_disorder.nf` — FASTA-последовательность и disorder-процент для одного белка, удобен для проверки.

## Установка Nextflow

```bash
conda create --name nf-env bioconda::nextflow
source activate nf-env
nextflow info


```

## RCSB PDBx/mmCIF

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
  rcsb/
    ab/
      1abc.cif.gz
    zz/
      9zzz.cif.gz
```

Например, структура `1abc` будет сохранена как:

```text
data/raw/rcsb/ab/1abc.cif.gz
```

## DisProt

Скачать текущую выгрузку DisProt в TSV:

```bash
nextflow run pipeline/download_disprot.nf
```

По умолчанию используется ссылка:

```text
https://disprot.org/api/v2/download?format=tsv&release=current&term_ontology=IDPO&term_ontology=GO
```

Файл будет сохранен сюда:

```text
data/raw/disprot/disprot_current_idpo_go.tsv
```

Переопределить URL или имя файла:

```bash
nextflow run pipeline/download_disprot.nf \
  --disprot_url 'https://disprot.org/api/v2/download?format=tsv&release=current&term_ontology=IDPO&term_ontology=GO' \
  --disprot_filename disprot_current_idpo_go.tsv
```

## UniProt последовательность + disorder процент

Собрать датасет для всех уникальных белков из DisProt:

```bash
nextflow run pipeline/build_disprot_sequence_disorder_dataset.nf
```

Для быстрой проверки можно ограничить число белков:

```bash
nextflow run pipeline/build_disprot_sequence_disorder_dataset.nf --limit 100
```

Пайплайн читает DisProt TSV отсюда:

```text
data/raw/disprot/disprot_current_idpo_go.tsv
```

Скачанные FASTA сохраняются сюда:

```text
data/raw/uniprot/
```

Итоговый датасет будет здесь:

```text
data/processed/disprot/sequence_disorder_dataset.tsv
```

В итоговой таблице есть:

- `UniProt ACC` — идентификатор белка.
- `Sequence` — полная аминокислотная последовательность из UniProt.
- `Protein Disorder Content` — доля disorder из DisProt.
- `Disorder percent` — доля disorder в процентах.
- `Disorder structural regions` — количество disorder-регионов по DisProt.

Для проверки одного конкретного белка можно использовать отдельный пайплайн:

```bash
nextflow run pipeline/download_uniprot_sequence_with_disorder.nf --uniprot_acc P03265
```

По умолчанию пайплайн читает DisProt TSV отсюда:

```text
data/raw/disprot/disprot_current_idpo_go.tsv
```

Результаты будут сохранены сюда:

```text
data/raw/uniprot/P03265.fasta
data/processed/disprot/P03265_disorder_summary.tsv
```

В summary-файле есть:

- `Sequence length FASTA` — длина полной последовательности из UniProt.
- `Protein Disorder Content` — доля disorder из DisProt.
- `Disorder percent` — та же доля в процентах.
- `Disorder structural regions` — сколько строк DisProt для этого белка имеют `Term namespace = Structural state` и `Term name = disorder`.

Для другого белка достаточно заменить UniProt accession:

```bash
nextflow run pipeline/download_uniprot_sequence_with_disorder.nf --uniprot_acc P49913
```
