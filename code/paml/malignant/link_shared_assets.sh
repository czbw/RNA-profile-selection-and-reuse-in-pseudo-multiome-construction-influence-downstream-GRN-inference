#!/usr/bin/env bash
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED="$AMLGRN_PROJECT/paml/healthy"      # = /path/to/AML-GRN

echo ">> bundle : $ROOT"
echo ">> shared : $SHARED"

# 1) data/ must be a symlink to the shared data tree (fragments, atac_meta, tar and rna come with
#    it). A stray real data/ directory would shadow the shared tree and hide the staged RNA
#    h5ads; detect and fix that without deleting data.
[ -d "$SHARED/data" ] || { echo "ERROR: $SHARED/data not found (Stage-1 data tree)."; exit 1; }
if [ -L "$ROOT/data" ]; then
  tgt="$(readlink -f "$ROOT/data")"
  if [ "$tgt" = "$(readlink -f "$SHARED/data")" ]; then
    echo "  data -> $tgt (symlink OK)"
  else
    echo "  [WARN] data symlink points at $tgt, not $SHARED/data — repointing"; rm -f "$ROOT/data"; ln -s "$SHARED/data" "$ROOT/data"; echo "  repointed data -> $SHARED/data"
  fi
elif [ -e "$ROOT/data" ]; then
  # a real path exists where the symlink should be.
  if [ -z "$(find "$ROOT/data" -mindepth 1 -print -quit 2>/dev/null)" ]; then
    rmdir "$ROOT/data" 2>/dev/null && ln -s "$SHARED/data" "$ROOT/data" && echo "  replaced empty real data/ with symlink -> $SHARED/data"
  else
    echo "  [ACTION NEEDED] $ROOT/data is a NON-EMPTY real directory shadowing $SHARED/data."
    echo "                  Its contents:"; find "$ROOT/data" -maxdepth 2 -mindepth 1 | sed 's/^/                    /'
    echo "                  If those files also live in $SHARED/data, remove the stray dir and symlink:"
    echo "                      rm -rf \"$ROOT/data\" && ln -s \"$SHARED/data\" \"$ROOT/data\""
    echo "                  If any RNA h5ad exists ONLY here, move it first:"
    echo "                      mv \"$ROOT/data/rna/\"malignant_*.h5ad \"$SHARED/data/rna/\" && rm -rf \"$ROOT/data\" && ln -s \"$SHARED/data\" \"$ROOT/data\""
  fi
else
  ln -s "$SHARED/data" "$ROOT/data"; echo "  linked data -> $SHARED/data"
fi

# 2) link the two large cisTarget feathers into the bundle's real databases/ (small sidecars stay in place).
mkdir -p "$ROOT/databases"
missing_db=0
for f in hg38_screen_v10_clust.regions_vs_motifs.rankings.feather \
         hg38_screen_v10_clust.regions_vs_motifs.scores.feather; do
  if [ -e "$ROOT/databases/$f" ]; then echo "  db feather present: $f"; continue; fi   # already there/linked — fine
  src=""
  for cand in "$SHARED/databases/$f" "$SHARED/$f"; do [ -f "$cand" ] && { src="$cand"; break; }; done
  if [ -n "$src" ]; then ln -s "$src" "$ROOT/databases/$f"; echo "  linked databases/$f -> $src"
  else echo "  [WARN] feather not found in $SHARED/databases or $SHARED: $f"; missing_db=1; fi
done
[ "$missing_db" -eq 0 ] || echo "  [note] a feather is missing — check the Stage-1 databases/ location"

# 3) check that the RNA h5ads are staged, at the shared location (independent of the bundle symlink).
for s in diagnosis relapse; do
  f="$SHARED/data/rna/malignant_${s}.h5ad"
  [ -f "$f" ] && echo "  RNA ok: malignant_${s}.h5ad ($(du -h "$f" | cut -f1))" \
              || echo "  [TODO] stage malignant_${s}.h5ad into $SHARED/data/rna (see MANIFEST)"
done

echo ">> done. Next: sbatch submit_bigmem.sh   (after the RNA h5ads are staged)"
