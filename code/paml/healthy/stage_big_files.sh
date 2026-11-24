#!/usr/bin/env bash
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${AMLGRN_SRC:-$AMLGRN_PROJECT}"          # source repo (override with AMLGRN_SRC if elsewhere)

echo ">> staging RNA h5ads -> data/rna/"
cp -v "$SRC"/data/paml/healthy_diagnosis.h5ad "$HERE/data/rna/"
cp -v "$SRC"/data/paml/healthy_remission.h5ad "$HERE/data/rna/"
cp -v "$SRC"/data/paml/healthy_relapse.h5ad   "$HERE/data/rna/"

echo ">> staging cisTarget feathers (46 GB) -> databases/"
cp -v "$SRC"/data/reference/motifs/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather "$HERE/databases/"
cp -v "$SRC"/data/reference/motifs/hg38_screen_v10_clust.regions_vs_motifs.scores.feather   "$HERE/databases/"

echo ">> staging ATAC tar (49 GB) -> data/"
cp -v "$SRC"/GSE235308_RAW.tar "$HERE/data/"

echo ">> staged. bundle is now self-contained. Verify with: du -sh $HERE"
echo ">> upload:  rsync -av --progress $HERE  <user>@cubic:<dest>/"
