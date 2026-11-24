#!/usr/bin/env bash
set -uo pipefail

ARM="${1:?usage: preflight_arm.sh <ARM> [STAGE] [DATASET] [SEED]}"
STAGE="${2:-pbmc10k}"; DS="${3:-10k}"; SEED="${4:-555}"
BUNDLE="${AMLGRN_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
OUT="$BUNDLE/out"; SCR="$BUNDLE/scripts"
FAIL=0
ok()   { printf '  [ OK ] %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*" >&2; FAIL=$((FAIL+1)); }

echo "== preflight: arm='$ARM' stage='$STAGE' dataset='$DS' seed='$SEED'"
echo "== bundle : $BUNDLE"

# ---- 1. the scripts exist -------------------------------------------------
for f in "$SCR/run_grn.sh" \
         "$SCR/chain_guard.sh"; do
  [ -s "$f" ] && ok "present: ${f#$BUNDLE/}" || bad "MISSING: ${f#$BUNDLE/}"
done

RESOLVE_OUT="$(AMLGRN_RESOLVE_ONLY=1 AMLGRN_ROOT="$BUNDLE" \
                bash "$SCR/run_grn.sh" "$STAGE" "$ARM" 2>&1)"
if [ $? -eq 0 ] && [ "${RESOLVE_OUT#resolve:}" != "$RESOLVE_OUT" ]; then
  ok "run_grn.sh ${RESOLVE_OUT}"
elif printf '%s' "$RESOLVE_OUT" | grep -q "unknown pairing"; then
  bad "run_grn.sh does NOT resolve '$ARM' -- it would die 'unknown pairing $ARM'. Upload the current scripts/run_grn.sh."
else
  bad "run_grn.sh resolve check failed for '$ARM': ${RESOLVE_OUT}"
fi

IS_METACELL=0
case "$ARM" in pseudomultiome|pseudomultiome_k*) IS_METACELL=1 ;; esac

# Assignment tables are needed only for single-cell constructions.
if [ "$IS_METACELL" -eq 0 ]; then
  CSV="$BUNDLE/uploads/arm_${ARM}_${DS}_seed${SEED}_full.csv"
  [ -s "$CSV" ] && ok "assignment table present: ${CSV#$BUNDLE/}" \
                || bad "MISSING assignment table $CSV"
fi

# ---- 4. inputs the arm will actually consume ------------------------------
if [ "$IS_METACELL" -eq 1 ]; then
  K="${ARM#pseudomultiome_k}"; [ "$K" = "$ARM" ] && K=5
  G="$OUT/native_rna_${STAGE}.metacell.h5ad"
  [ -s "$G" ] && ok "metacell GEX present: ${G#$BUNDLE/}" || bad "MISSING $G"
  [ -s "$SCR/metacell_floor.py" ] \
    && ok "metacell_floor.py present (k=$K floor will be asserted at runtime)" \
    || bad "MISSING scripts/metacell_floor.py -- run_grn.sh calls it and will die"
fi
for f in "$OUT/cistopic_obj_${STAGE}.pkl" "$OUT/menr_cistarget_${STAGE}.hdf5" \
         "$OUT/atac_chain_${STAGE}.sha256"; do
  [ -s "$f" ] && ok "chain artifact: ${f#$BUNDLE/}" || bad "MISSING chain artifact $f"
done
[ -d "$OUT/region_sets_${STAGE}" ] && ok "region_sets_${STAGE}/" || bad "MISSING $OUT/region_sets_${STAGE}"

# ---- 5. would this arm overwrite a completed one? --------------------------
if [ -s "$OUT/eRegulon_AUC_${STAGE}_${ARM}.h5mu" ]; then
  ok "NOTE: '$ARM' already complete -- run_grn.sh will skip it (delete outputs to force a rerun)"
fi

echo
if [ "$FAIL" -ne 0 ]; then
  echo "== PREFLIGHT FAILED ($FAIL problem(s)) -- do NOT sbatch until these are fixed." >&2
  exit 1
fi
echo "== PREFLIGHT PASSED. Run scripts/run_grn.sh $STAGE $ARM in the configured environment with AMLGRN_ROOT=$BUNDLE."
