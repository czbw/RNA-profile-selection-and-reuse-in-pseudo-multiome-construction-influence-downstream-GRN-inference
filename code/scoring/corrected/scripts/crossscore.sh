#!/usr/bin/env bash
set -uo pipefail

S2="${AMLGRN_S2:?export AMLGRN_S2}"          # paml/malignant/out
OUT="${AMLGRN_S3_OUT:?export AMLGRN_S3_OUT}" # scoring/corrected/out
ARM="${AMLGRN_S3_PAIRING:?export AMLGRN_S3_PAIRING}"
mkdir -p "$OUT"

partof() { local f="$1"; printf '%s.part.%s' "${f%.*}" "${f##*.}"; }

run_one () {
  local CELLS="$1" DEFS="$2"
  local MUD="$S2/mudata_${CELLS}_${ARM}.h5mu"
  local EREG="$S2/eRegulons_${DEFS}_${ARM}.tsv"
  local AUC="$OUT/eRegulon_AUC_${CELLS}_${ARM}_def${DEFS}.h5mu"

  echo "===== [${ARM}: ${CELLS} cells x ${DEFS} defs] ====="
  [ -s "$MUD" ]  || { echo "FAILED: missing $MUD";  return 1; }
  [ -s "$EREG" ] || { echo "FAILED: missing $EREG"; return 1; }
  if [ -s "$AUC" ]; then echo "  already present — skip"; return 0; fi

  scenicplus grn_inference AUCell \
      --eRegulon_fname "$EREG" \
      --multiome_mudata_fname "$MUD" \
      --aucell_out_fname "$(partof "$AUC")" \
      --n_cpu 1 || { echo "FAILED: AUCell ${ARM} ${CELLS}x${DEFS}"; return 1; }

  mv -f "$(partof "$AUC")" "$AUC" || { echo "FAILED: commit ${ARM} ${CELLS}x${DEFS}"; return 1; }
  echo "  -> $AUC"
  return 0
}

FAIL=0
# Diagnosis cells are the memory peak: mudata_Diagnosis_repair carries 112,556 cells against
# Relapse's 70,876, over the same 624,963-region dense + rank matrix. Run Diagnosis first so a
# mis-sized job fails early rather than after the cheaper half.
run_one Diagnosis Relapse   || FAIL=$((FAIL+1))
run_one Relapse   Diagnosis || FAIL=$((FAIL+1))

echo
echo "==== STAGE 3 CROSS-SCORING SUMMARY (${ARM}) ===="
for f in "$OUT/eRegulon_AUC_Diagnosis_${ARM}_defRelapse.h5mu" \
         "$OUT/eRegulon_AUC_Relapse_${ARM}_defDiagnosis.h5mu"; do
  if [ -s "$f" ]; then printf '  OK      %s  (%s)\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else                 printf '  MISSING %s\n' "$(basename "$f")"; fi
done
ls "$OUT"/*.part.* >/dev/null 2>&1 && echo "  NOTE: stray .part.* present (a killed run; safe to delete)"

[ "$FAIL" -eq 0 ] || { echo "==== $FAIL run(s) FAILED ===="; exit 1; }
echo "==== CROSS-SCORING COMPLETE — Gene_based is the only valid modality here ===="
