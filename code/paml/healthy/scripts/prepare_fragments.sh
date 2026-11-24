#!/usr/bin/env bash
set -uo pipefail
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT}"
TAR="$ROOT/data/GSE235308_RAW.tar"
RAW="$ROOT/data/fragments_raw"
OUT="$ROOT/data/fragments"
TMP="$ROOT/work/.frag_tmp"
mkdir -p "$RAW" "$OUT" "$TMP"
command -v bgzip >/dev/null || { echo "ERROR: bgzip not found (conda activate scenicplus)"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found"; exit 1; }

# (1) extract filtered fragments from the tar if not already extracted
if [ -f "$TAR" ]; then
  echo ">> extracting *_filtered_fragments.tsv.gz from $(basename "$TAR") -> $RAW"
  # -C must precede the member pattern: a non-permuting tar (POSIXLY_CORRECT GNU tar,
  # or bsdtar in minimal images) otherwise treats '-C $RAW' as extra member patterns and
  # extracts into the CWD instead of $RAW. Try GNU-tar syntax first, fall back to bsdtar.
  tar -x -v -f "$TAR" -C "$RAW" --wildcards --no-anchored '*_filtered_fragments.tsv.gz' 2>/dev/null \
    || tar -x -v -f "$TAR" -C "$RAW" '*_filtered_fragments.tsv.gz' \
    || true
else
  echo ">> NOTE: $TAR not present. Assuming $RAW already holds *_filtered_fragments.tsv.gz."
fi

shopt -s nullglob
frags=( "$RAW"/*_filtered_fragments.tsv.gz )
[[ ${#frags[@]} -gt 0 ]] || { echo "ERROR: no *_filtered_fragments.tsv.gz in $RAW"; exit 1; }
echo ">> ${#frags[@]} samples to convert"

for fp in "${frags[@]}"; do
  base="$(basename "$fp" _filtered_fragments.tsv.gz)"     # e.g. GSM7498708_AML16_DX
  out="${OUT}/${base}.fragments.tsv.gz"
  if [[ -f "${out}.tbi" ]]; then echo ">> $base already done"; continue; fi
  echo ">> $base"
  zcat "$fp" \
    | awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1, $2-1, $3, $6, 1}' \
    | sort -k1,1 -k2,2n -S 4G -T "$TMP" \
    | bgzip -c > "$out"
  tabix -p bed -f "$out"
done
rm -rf "$TMP"
echo ">> done. converted fragments + .tbi in $OUT"
