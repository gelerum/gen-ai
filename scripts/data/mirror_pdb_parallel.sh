#!/bin/bash
# Параллельное зеркалирование mmCIF-архива PDB в data/raw/.
# rsync однопоточный, поэтому распараллеливаем по подкаталогам (hh/, bn/, ...).
#
#   scripts/data/mirror_pdb_parallel.sh [DIR] [JOBS]
#
# DIR  — куда складывать (по умолчанию <repo>/data/raw/pdb_mmCIF,
#        независимо от текущего рабочего каталога).
# JOBS — число параллельных потоков rsync (по умолчанию 6, больше 8 не ставь).

set -euo pipefail

# Корень репозитория: этот скрипт лежит в <repo>/scripts/data/.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

MIRRORDIR="${1:-${REPO_ROOT}/data/raw/pdb_mmCIF}"
JOBS="${2:-6}"                   # больше 8 не ставь: сервер общий
SERVER="rsync.rcsb.org"
PORT=33444
SRC="${SERVER}::ftp_data/structures/divided/mmCIF"

mkdir -p "$MIRRORDIR"

echo "Зеркалирую в: $MIRRORDIR"
echo "Получаю список подкаталогов..."
mapfile -t DIRS < <(
  rsync --list-only --port="$PORT" "${SRC}/" \
    | awk '$1 ~ /^d/ && $NF != "." { print $NF }'
)
echo "Подкаталогов: ${#DIRS[@]}, потоков: $JOBS"

sync_one() {
  local d="$1"
  mkdir -p "${MIRRORDIR}/${d}"
  # --delete здесь безопасен: скоуп ограничен одним подкаталогом
  rsync -rlpt --delete --port="$PORT" \
    "${SRC}/${d}/" "${MIRRORDIR}/${d}/"
}
export -f sync_one
export MIRRORDIR SRC PORT

printf '%s\n' "${DIRS[@]}" \
  | xargs -P "$JOBS" -I{} bash -c 'sync_one "$@"' _ {}

echo
echo "Файлов: $(find "$MIRRORDIR" -name '*.cif.gz' | wc -l)"
echo "Объём:  $(du -sh "$MIRRORDIR" | cut -f1)"
