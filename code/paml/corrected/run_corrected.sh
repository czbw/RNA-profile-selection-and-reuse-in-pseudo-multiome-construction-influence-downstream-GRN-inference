#!/usr/bin/env bash
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus

ROOT="${AMLGRN_ROOT:?}"                       # = paml/malignant (chain + out/ + work/)
BUNDLE="${AMLGRN_REPAIR_BUNDLE:?}"            # = paml/corrected
STAGES="${AMLGRN_STAGES:-Diagnosis Relapse}"
OUT="$ROOT/out"; WORK="$ROOT/work"

echo "=================================================================="
echo " AML RE-PAIRING   $(date)"
echo " python  : $(which python)"
echo " chain   : $OUT"
echo " stages  : $STAGES   ncpu: ${AMLGRN_NCPU:-20}"
echo " TMPDIR  : ${TMPDIR:-unset}"
echo "=================================================================="

FAIL=0

# ---- [1] assert the ATAC chain BEFORE anything is written -----------------------------
echo; echo "---- [1/4] verifying the shared ATAC chain (never rebuilt) ----"
for STAGE in $STAGES; do
  for f in "$OUT/cistopic_obj_${STAGE}.pkl" "$OUT/menr_cistarget_${STAGE}.hdf5"; do
    [ -s "$f" ] || { echo "  MISSING CHAIN ARTIFACT $f -- refusing to rebuild it here"; FAIL=1; }
  done
  [ -d "$OUT/region_sets_${STAGE}" ] || { echo "  MISSING $OUT/region_sets_${STAGE}"; FAIL=1; }
done
[ "$FAIL" -eq 0 ] || { echo "==== ABORTING: chain incomplete ===="; exit 1; }
echo "  chain present for: $STAGES"

# ---- [2] stage EVERY table before the builder runs ------------------------------------
echo; echo "---- [2/4] staging repaired cell_data tables ----"
for STAGE in $STAGES; do
  SRC="$BUNDLE/uploads/cell_data_${STAGE}.repaired.csv"
  [ -s "$SRC" ] || { echo "  MISSING $SRC (scp it from paml/corrected/out/)"; FAIL=1; continue; }
  cp -f "$SRC" "$WORK/cell_data_${STAGE}.repaired.csv" || { echo "  cp FAILED for $STAGE"; FAIL=1; }
  echo "  $STAGE: $(wc -l < "$SRC") lines -> work/cell_data_${STAGE}.repaired.csv"
done
[ "$FAIL" -eq 0 ] || { echo "==== ABORTING: could not stage all tables ===="; exit 1; }

# ---- [3] build the paired RNA objects, ONCE for all staged stages ---------------------
echo; echo "---- [3/4] build_paired_rna.py (all staged stages in one call) ----"
NEED_BUILD=0
for STAGE in $STAGES; do
  [ -s "$OUT/paired_rna_${STAGE}.repair.h5ad" ] || NEED_BUILD=1
done
if [ "$NEED_BUILD" -eq 0 ]; then
  echo "  all paired_rna_*.repair.h5ad present -- skip"
else
  ( cd "$BUNDLE/scripts" && python build_paired_rna.py ) \
    || { echo "  build_paired_rna_repair FAILED"; echo "==== ABORTING ===="; exit 1; }
  for STAGE in $STAGES; do
    [ -s "$OUT/paired_rna_${STAGE}.repair.h5ad" ] \
      || { echo "  builder finished but paired_rna_${STAGE}.repair.h5ad is missing/empty"; FAIL=1; }
  done
  [ "$FAIL" -eq 0 ] || { echo "==== ABORTING: builder did not produce every stage ===="; exit 1; }
fi

# ---- [4] GRN half per stage, under TAG <STAGE>_repair ---------------------------------
echo; echo "---- [4/4] GRN half ----"
for STAGE in $STAGES; do
  if [ -s "$OUT/eRegulon_AUC_${STAGE}_repair.h5mu" ]; then
    echo "  [$STAGE] already complete -- skip"
    continue
  fi
  t0=$(date +%s)
  ( cd "$BUNDLE/scripts" && bash run_grn.sh "$STAGE" repair ) \
    || { echo "  [$STAGE] GRN FAILED"; FAIL=$((FAIL+1)); continue; }
  echo "  [$STAGE] GRN half took $(( ($(date +%s) - t0) / 60 )) min"
done

echo
echo "==== RE-PAIRING SUMMARY  $(date) ===="
for STAGE in $STAGES; do
  wc -l "$OUT/eRegulons_${STAGE}_repair.tsv" 2>/dev/null || echo "  (no eRegulons for $STAGE)"
done
if [ "$FAIL" -ne 0 ]; then echo "  $FAIL step(s) FAILED"; echo "==== INCOMPLETE ===="; exit 1; fi
echo "==== RE-PAIRING COMPLETE ===="
