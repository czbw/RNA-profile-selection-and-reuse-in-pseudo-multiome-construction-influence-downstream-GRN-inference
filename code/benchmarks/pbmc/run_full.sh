#!/usr/bin/env bash
set -uo pipefail

BUNDLE="${AMLGRN_ROOT:?submit_full.sh must export AMLGRN_ROOT}"
ROOT="$BUNDLE"
DS="${1:-10k}"
CFG=config.yaml
OUT="$BUNDLE/out"
WORK="$BUNDLE/work"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus
export AMLGRN_MACS="$(conda run -n macs3env which macs3)"
export AMLGRN_CONFIG="$CFG"
export AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-20}"

cd "$BUNDLE/scripts" || { echo "FATAL: no $BUNDLE/scripts"; exit 3; }
STAGE=$(python -c "import yaml;print(yaml.safe_load(open('$CFG'))['stages'][0])")
export AMLGRN_STAGE="$STAGE"

echo ">> python  : $(which python)"
echo ">> macs3   : $AMLGRN_MACS"
echo ">> cwd     : $(pwd)   (bare-filename convention)"
echo ">> stage   : $STAGE   dataset: $DS   ncpu: $AMLGRN_NCPU"
echo ">> DB      : ${AMLGRN_DB:-$BUNDLE/databases}"
echo ">> pairings: truth link pseudomultiome"

FAIL=0
step () {
  local n="$1" label="$2"; shift 2
  echo; echo "---- [$n] $label ----"
  if "$@"; then echo "---- [$n] ok"; else echo "---- [$n] FAILED"; FAIL=$((FAIL+1)); fi
}

# ---- prepare / annotate / assemble the three arms --------------------------
[ -s "$WORK/$DS/gex.h5ad" ] \
  && echo "[1] gex.h5ad present -- skip" \
  || step 1 "prepare"  python prepare.py --config "$CFG" --dataset "$DS"

[ "$FAIL" -eq 0 ] && { [ -s "$WORK/$DS/gex_annotated.h5ad" ] \
  && echo "[2] gex_annotated.h5ad present -- skip" \
  || step 2 "annotate" python annotate.py --config "$CFG" --dataset "$DS" --resolution 1.0; }

[ "$FAIL" -eq 0 ] && { [ -s "$OUT/paired_rna_${STAGE}.link.h5ad" ] \
  && echo "[3] pairings present -- skip" \
  || step 3 "make pairings (truth / link / pseudomultiome)" \
       python make_pairings.py --config "$CFG" --dataset "$DS"; }

[ "$FAIL" -eq 0 ] && { [ -s "$OUT/cistopic_obj_${STAGE}.pkl" ] \
  && echo "[4] cisTopic object present -- skip" \
  || step 4 "ATAC chain (pseudobulk -> peaks -> cisTopic -> topics)" python atac_chain.py; }

[ "$FAIL" -eq 0 ] && { [ -d "$OUT/region_sets_${STAGE}" ] \
  && echo "[5] region sets present -- skip" \
  || step 5 "region sets" python region_sets.py; }

[ "$FAIL" -eq 0 ] && { [ -s "$OUT/menr_cistarget_${STAGE}.hdf5" ] \
  && echo "[5b] motif enrichment present -- skip" \
  || step 5b "motif enrichment (cisTarget, ONCE)" \
       bash -c 'scenicplus grn_inference motif_enrichment_cistarget \
         --region_set_folder "'"$OUT"'/region_sets_'"$STAGE"'" \
         --cistarget_db_fname "'"${AMLGRN_DB:-$ROOT/databases}"'/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather" \
         --output_fname_cistarget_result "'"$OUT"'/menr_cistarget_'"$STAGE"'.hdf5" \
         --temp_dir "'"$ROOT"'/work/grn_tmp_'"$STAGE"'_shared" --species homo_sapiens \
         --path_to_motif_annotations "'"${AMLGRN_DB:-$ROOT/databases}"'/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
         --annotation_version v10nr_clust --n_cpu 1'; }

# Record ATAC-chain hashes so that any difference between arms can be attributed to the pairing
# and not to the ATAC side.
[ "$FAIL" -eq 0 ] && step 5c "record ATAC-chain hashes" \
  bash chain_guard.sh record "$STAGE"

# ---- GRN runs against that one frozen ATAC chain ---------------------------
# run_grn.sh resolves the GEX per pairing:
#   truth / link   -> out/paired_rna_<stage>.<arm>.h5ad   (multiome mode)
#   pseudomultiome -> out/native_rna_<stage>.metacell.h5ad
# Each arm re-verifies the chain at entry and aborts non-zero on any mismatch.
ARMS="${AMLGRN_ARMS:-truth link pseudomultiome}"
if [ "$FAIL" -eq 0 ]; then
  for ARM in $ARMS; do
    if [ -s "$OUT/eRegulon_AUC_${STAGE}_${ARM}.h5mu" ]; then
      echo "[6:$ARM] already complete -- skip"
    else
      step "6:$ARM" "GRN + AUCell ($ARM)" bash run_grn.sh "$STAGE" "$ARM"
    fi
  done
fi

# ---- [7] re-verify the chain after the last arm ----------------------------
[ "$FAIL" -eq 0 ] && step 7 "re-verify ATAC chain (post-arms)" \
  bash chain_guard.sh verify "$STAGE" "post-all-arms"

echo
echo "==== PBMC 10k SUMMARY ===="
if [ "$FAIL" -ne 0 ]; then
  echo "  $FAIL step(s) FAILED -- see above."
  echo "==== INCOMPLETE ===="; exit 1
fi
for ARM in truth link pseudomultiome; do
  f="$OUT/eRegulon_AUC_${STAGE}_${ARM}.h5mu"
  if [ -s "$f" ]; then printf '  OK      %s (%s)\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else printf '  MISSING %s\n' "$(basename "$f")"; FAIL=$((FAIL+1)); fi
done
[ "$FAIL" -eq 0 ] || { echo "==== INCOMPLETE -- expected outputs missing ===="; exit 1; }
echo
echo "  All three arms complete. Next: score link and pseudomultiome AGAINST truth"
echo "  (edge recovery, activity correlation, RNA-side collapse)."
echo "==== PBMC 10k COMPLETE ===="
