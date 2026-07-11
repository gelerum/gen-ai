# Загрузка RCSB PDB в формате PDBx/mmCIF

Проект скачивает структуры RCSB/wwPDB в формате PDBx/mmCIF (`.cif.gz`) и складывает их в `data/raw/rcsb`.

Пайплайн находится в `pipeline/download_pdb_mmcif.nf` и использует `curl` внутри задач Nextflow.

## Установка Nextflow

```bash
conda create --name nf-env bioconda::nextflow
source activate nf-env
nextflow info


```

## Запуск

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
