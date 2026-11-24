#!/usr/bin/env bash
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

DS="${1:-3k}"
ROOT="${AMLGRN_PROJECT}"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"

echo "=============================================================="
echo " PBMC link arm (ArchR)   dataset=$DS   $(date)"
echo "=============================================================="

[ -x "$MM" ] || { echo "FATAL: micromamba not at $MM" >&2; exit 3; }
"$MM" env list | grep -q '^\s*archr' || { echo "FATAL: micromamba env 'archr' not found" >&2; exit 3; }

for f in "$ROOT/data/pbmc/pbmc_granulocyte_sorted_${DS}_atac_fragments.tsv.gz" \
         "$ROOT/data/pbmc/pbmc_granulocyte_sorted_${DS}_filtered_feature_bc_matrix.h5" \
         "$ROOT/benchmarks/pbmc/scripts/archr_link.R"; do
  [ -s "$f" ] || { echo "FATAL: missing/empty $f" >&2; exit 3; }
done

# fail early if the R toolchain is incomplete, before starting the ArchR run
"$MM" run -n archr Rscript -e '
  need <- c("ArchR","Seurat","chromVARmotifs","rhdf5","jsonlite",
            "BSgenome.Hsapiens.UCSC.hg38")
  miss <- need[!sapply(need, requireNamespace, quietly=TRUE)]
  if (length(miss)) { cat("FATAL missing R packages:", paste(miss, collapse=", "), "\n"); quit(status=3) }
  cat(">> R", paste(R.version$major, R.version$minor, sep="."),
      "| ArchR", as.character(packageVersion("ArchR")),
      "| Seurat", as.character(packageVersion("Seurat")), "\n")
' || { echo "FATAL: R toolchain incomplete -- see above" >&2; exit 3; }

echo ">> launching archr_link.R"
"$MM" run -n archr Rscript "$ROOT/benchmarks/pbmc/scripts/archr_link.R" \
      --dataset "$DS" --root "$ROOT" --threads 6
rc=$?

OUT="$ROOT/benchmarks/pbmc/work/$DS/archr"
echo
if [ "$rc" -ne 0 ]; then
  echo "==== LINK ARM FAILED (rc=$rc) ===="
  echo "  If the failure is inside addGeneIntegrationMatrix on Seurat 5, the fix is a pinned"
  echo "  r-seurat=4.* env -- NOT a reimplementation of the link. Comparability with the AML"
  echo "  arm is the entire point."
  exit "$rc"
fi
for f in "$OUT/link_${DS}.csv" "$OUT/link_report_${DS}.json"; do
  if [ -s "$f" ]; then printf '  OK      %s (%s)\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else printf '  MISSING %s\n' "$f"; rc=1; fi
done
[ "$rc" -eq 0 ] || { echo "==== LINK ARM INCOMPLETE ===="; exit 1; }
echo "==== LINK ARM COMPLETE ===="
