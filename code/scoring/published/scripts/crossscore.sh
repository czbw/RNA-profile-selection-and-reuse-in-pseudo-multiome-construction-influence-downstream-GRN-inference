#!/usr/bin/env bash
set -uo pipefail

S2="${AMLGRN_S2:?export AMLGRN_S2}"          # paml/malignant/out
OUT="${AMLGRN_S3_OUT:?export AMLGRN_S3_OUT}" # scoring/published/out
mkdir -p "$OUT"

partof() { local f="$1"; printf '%s.part.%s' "${f%.*}" "${f##*.}"; }

run_one () {
  local CELLS="$1" DEFS="$2"
  local MUD="$S2/mudata_${CELLS}_link.h5mu"
  local EREG="$S2/eRegulons_${DEFS}_link.tsv"
  local AUC="$OUT/eRegulon_AUC_${CELLS}_link_def${DEFS}.h5mu"

  echo "===== [${CELLS} cells x ${DEFS} defs] ====="
  [ -s "$MUD" ]  || { echo "FAILED: missing $MUD";  return 1; }
  [ -s "$EREG" ] || { echo "FAILED: missing $EREG"; return 1; }
  if [ -s "$AUC" ]; then echo "  already present — skip"; return 0; fi

  scenicplus grn_inference AUCell \
      --eRegulon_fname "$EREG" \
      --multiome_mudata_fname "$MUD" \
      --aucell_out_fname "$(partof "$AUC")" \
      --n_cpu 1 || { echo "FAILED: AUCell ${CELLS}x${DEFS}"; return 1; }

  mv -f "$(partof "$AUC")" "$AUC" || { echo "FAILED: commit ${CELLS}x${DEFS}"; return 1; }
  echo "  -> $AUC"
  return 0
}

FAIL=0
# Diagnosis cells are the memory peak (100,659 cells x 624,963 regions dense + rank matrix);
# run them first so a mis-sized job fails early rather than after the cheap half.
run_one Diagnosis Relapse   || FAIL=$((FAIL+1))
run_one Relapse   Diagnosis || FAIL=$((FAIL+1))

echo
echo "==== STAGE 3 CROSS-SCORING SUMMARY ===="
for f in "$OUT/eRegulon_AUC_Diagnosis_link_defRelapse.h5mu" \
         "$OUT/eRegulon_AUC_Relapse_link_defDiagnosis.h5mu"; do
  if [ -s "$f" ]; then printf '  OK      %s  (%s)\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else                 printf '  MISSING %s\n' "$(basename "$f")"; fi
done
ls "$OUT"/*.part.* >/dev/null 2>&1 && echo "  NOTE: stray .part.* present (a killed run; safe to delete)"

[ "$FAIL" -eq 0 ] || { echo "==== $FAIL run(s) FAILED ===="; exit 1; }
echo "==== CROSS-SCORING COMPLETE — Gene_based is the only valid modality here ===="
