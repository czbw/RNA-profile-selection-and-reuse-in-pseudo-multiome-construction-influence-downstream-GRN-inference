#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW="${PROJECT_ROOT}/data/paml/fragments_raw"
OUT="${PROJECT_ROOT}/data/paml/fragments"
TMP="${PROJECT_ROOT}/preprocessing/work"
mkdir -p "$OUT" "$TMP"

command -v bgzip >/dev/null || { echo "ERROR: bgzip not found (conda activate scenicplus)"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found (conda activate scenicplus)"; exit 1; }

shopt -s nullglob
frags=( "$RAW"/*_filtered_fragments.tsv.gz )
[[ ${#frags[@]} -gt 0 ]] || { echo "ERROR: no fragments in $RAW — run extract_atac.sh"; exit 1; }

for fp in "${frags[@]}"; do
  base="$(basename "$fp" _filtered_fragments.tsv.gz)"     # e.g. GSM7498708_AML16_DX
  out="${OUT}/${base}.fragments.tsv.gz"
  if [[ -f "${out}.tbi" ]]; then echo ">> $base already done"; continue; fi
  echo ">> $base"
  # header-skip + GRanges->BED + sort + bgzip
  zcat "$fp" \
    | awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1, $2-1, $3, $6, 1}' \
    | sort -k1,1 -k2,2n -S 2G -T "$TMP" \
    | bgzip -c > "$out"
  tabix -p bed -f "$out"
done

rm -rf "$TMP"
echo ">> done. converted fragments + .tbi in $OUT"
echo ">> next: python make_cell_data.py"
