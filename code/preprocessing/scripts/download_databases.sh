#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DB="${PROJECT_ROOT}/data/reference/motifs"
mkdir -p "$DB"
cd "$DB"

echo ">> Genome build expected: hg38. Abort now (Ctrl-C) if the audit did NOT confirm hg38."

fetch() { echo ">> fetch $(basename "$1")"; wget -c "$1"; }

BASE="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based"
fetch "${BASE}/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
fetch "${BASE}/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
fetch "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
fetch "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/allTFs_hg38.txt"
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

# --- truncation sanity check (expected sizes, bytes) ---------------------------------------
declare -A MIN=(
  [hg38_screen_v10_clust.regions_vs_motifs.rankings.feather]=35000000000   # ~32.8 GB
  [hg38_screen_v10_clust.regions_vs_motifs.scores.feather]=13000000000     # ~12.9 GB
)
bad=0
for f in "${!MIN[@]}"; do
  have=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if (( have < MIN[$f] )); then
    echo "!! INCOMPLETE: $f is $have bytes (< ${MIN[$f]}). Re-run this script to resume."; bad=1
  fi
done
echo ">> databases in $DB"; ls -lh "$DB"
(( bad == 0 )) && echo ">> ALL DBs complete." || { echo ">> Some DBs incomplete — re-run to finish."; exit 1; }
