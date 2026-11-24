#!/usr/bin/env bash
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus

HOUT="${AMLGRN_HEALTHY_OUT:?}"        # CUBIC: $ROOT/out      (healthy chain + healthy outputs)
HWORK="${AMLGRN_HEALTHY_WORK:?}"      # CUBIC: $ROOT/work
HRNA="${AMLGRN_HEALTHY_RNA:?}"        # CUBIC: $ROOT/data/rna
MALIG="${AMLGRN_MALIG:?}"             # paml/malignant/out
BUNDLE="${AMLGRN_HR_BUNDLE:?}"        # paml/healthy_corrected
OUT="${AMLGRN_HR_OUT:?}"              # paml/healthy_corrected/out
STAGES="${AMLGRN_STAGES:-Diagnosis Relapse}"
DEFS="${AMLGRN_DEFS:-Diagnosis Relapse}"

echo "=================================================================="
echo " HEALTHY RE-PAIRING   $(date)"
echo " python      : $(which python)"
echo " healthy out : $HOUT"
echo " malignant   : $MALIG"
echo " stages      : $STAGES   defs: $DEFS   ncpu: ${AMLGRN_NCPU:-8}"
echo " TMPDIR      : ${TMPDIR:-unset}"
echo "=================================================================="

FAIL=0

# ---- [1] assert everything BEFORE anything is written ----------------------------------
echo; echo "---- [1/4] verifying the healthy ATAC chain and the malignant repair vocabularies ----"
for STAGE in $STAGES; do
  f="$HOUT/cistopic_obj_${STAGE}.pkl"
  if [ -s "$f" ]; then
    echo "  chain OK : $f"
  else
    echo "  MISSING CHAIN ARTIFACT $f"
    echo "           This is the fork in the road. If the healthy cisTopic object is genuinely"
    echo "           absent, it must be rebuilt by the Arm-A run at \$ROOT (atac_chain.py) and that is a"
    echo "           MULTI-DAY job, not this one. Refusing to rebuild it implicitly here."
    FAIL=1
  fi
done
for D in $DEFS; do
  f="$MALIG/eRegulons_${D}_repair.tsv"
  [ -s "$f" ] && echo "  defs  OK : $f" || { echo "  MISSING $f (the malignant repair vocabulary)"; FAIL=1; }
done
for STAGE in $STAGES; do
  f="$HRNA/healthy_$(echo "$STAGE" | tr '[:upper:]' '[:lower:]').h5ad"
  [ -s "$f" ] || { echo "  MISSING healthy RNA $f"; FAIL=1; }
done
# the builder searches all three healthy RNA files regardless of the stages requested
f="$HRNA/healthy_remission.h5ad"
[ -s "$f" ] || { echo "  MISSING healthy RNA $f (needed for the donor lookup even when Remission is not built)"; FAIL=1; }
[ "$FAIL" -eq 0 ] || { echo "==== ABORTING: prerequisites incomplete ===="; exit 1; }

# ---- [2] stage EVERY table before the builder runs -------------------------------------
echo; echo "---- [2/4] staging repaired healthy cell_data tables ----"
mkdir -p "$HWORK"
for STAGE in $STAGES; do
  SRC="$BUNDLE/uploads/cell_data_${STAGE}.healthy_repaired.csv"
  [ -s "$SRC" ] || { echo "  MISSING $SRC (scp it from paml/healthy_corrected/out/)"; FAIL=1; continue; }
  cp -f "$SRC" "$HWORK/cell_data_${STAGE}.healthy_repaired.csv" \
    || { echo "  cp FAILED for $STAGE"; FAIL=1; continue; }
  echo "  $STAGE: $(wc -l < "$SRC") lines -> work/cell_data_${STAGE}.healthy_repaired.csv"
done
[ "$FAIL" -eq 0 ] || { echo "==== ABORTING: could not stage all tables ===="; exit 1; }

# ---- [3] build the paired healthy GEX, ONCE for all staged stages ----------------------
echo; echo "---- [3/4] build_paired_rna.py (all staged stages in one call) ----"
NEED_BUILD=0
for STAGE in $STAGES; do
  [ -s "$HOUT/paired_rna_${STAGE}.healthy_repair.h5ad" ] || NEED_BUILD=1
done
if [ "$NEED_BUILD" -eq 0 ]; then
  echo "  all paired_rna_*.healthy_repair.h5ad present -- skip"
else
  ( cd "$BUNDLE/scripts" && python build_paired_rna.py ) \
    || { echo "  build_paired_rna_healthy_repair FAILED"; echo "==== ABORTING ===="; exit 1; }
  for STAGE in $STAGES; do
    [ -s "$HOUT/paired_rna_${STAGE}.healthy_repair.h5ad" ] \
      || { echo "  builder finished but paired_rna_${STAGE}.healthy_repair.h5ad is missing/empty"; FAIL=1; }
  done
  [ "$FAIL" -eq 0 ] || { echo "==== ABORTING: builder did not produce every stage ===="; exit 1; }
fi

# ---- [4] mudata + the cross-scored 2x2 -------------------------------------------------
echo; echo "---- [4/4] mudata rebuild + cross-scoring ----"
t0=$(date +%s)
( cd "$BUNDLE/scripts" && bash crossscore.sh ) || FAIL=$((FAIL+1))
echo "  cross-scoring took $(( ($(date +%s) - t0) / 60 )) min"

echo
echo "==== HEALTHY RE-PAIRING SUMMARY  $(date) ===="
for HCELLS in $STAGES; do
  for MDEFS in $DEFS; do
    f="$OUT/eRegulon_AUC_H${HCELLS}_repair_defM${MDEFS}.h5mu"
    if [ -s "$f" ]; then printf '  OK      %s  (%s)\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
    else                 printf '  MISSING %s\n' "$(basename "$f")"; fi
  done
done
if [ "$FAIL" -ne 0 ]; then echo "  $FAIL step(s) FAILED"; echo "==== INCOMPLETE ===="; exit 1; fi
echo "==== HEALTHY RE-PAIRING COMPLETE ===="
