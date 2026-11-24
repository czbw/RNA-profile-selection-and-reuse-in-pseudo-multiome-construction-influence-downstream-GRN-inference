#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"        # repo root (AML-GRN)
TAR="${ROOT}/GSE235308_RAW.tar"
RAW="${PROJECT_ROOT}/data/paml/fragments_raw"
META="${PROJECT_ROOT}/data/paml/atac_meta"
mkdir -p "$RAW" "$META"

[[ -f "$TAR" ]] || { echo "ERROR: $TAR not found"; exit 1; }

PATTERNS=( '*_filtered_fragments.tsv.gz' '*_filtered_barcodes.tsv.gz' )
if [[ "${1:-}" == "--with-peak-matrix" ]]; then
  PATTERNS+=( '*_filtered_peak_matrix.mtx.gz' )
fi

# metadata (idempotent — skip if already present)
if [[ "$(ls -1 "$META"/*_filtered_metadata.tsv.gz 2>/dev/null | wc -l)" -lt 64 ]]; then
  echo ">> extracting metadata -> $META"
  tar --wildcards -xvf "$TAR" -C "$META" '*_filtered_metadata.tsv.gz'
fi

echo ">> extracting fragments/barcodes -> $RAW  (this reads most of the 52 GB; be patient)"
tar --wildcards -xvf "$TAR" -C "$RAW" "${PATTERNS[@]}"

echo ">> extracted:"
echo "   fragments: $(ls -1 "$RAW"/*_filtered_fragments.tsv.gz 2>/dev/null | wc -l)"
echo "   barcodes : $(ls -1 "$RAW"/*_filtered_barcodes.tsv.gz 2>/dev/null | wc -l)"
echo "   metadata : $(ls -1 "$META"/*_filtered_metadata.tsv.gz 2>/dev/null | wc -l)"
echo ">> next: bash prepare_fragments.sh"
