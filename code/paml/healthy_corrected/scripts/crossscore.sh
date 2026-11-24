#!/usr/bin/env bash
set -uo pipefail

HOUT="${AMLGRN_HEALTHY_OUT:?export AMLGRN_HEALTHY_OUT}"   # CUBIC: $ROOT/out (the healthy-cell run owns the root)
MALIG="${AMLGRN_MALIG:?export AMLGRN_MALIG}"              # paml/malignant/out
OUT="${AMLGRN_HR_OUT:?export AMLGRN_HR_OUT}"              # paml/healthy_corrected/out
STAGES="${AMLGRN_STAGES:-Diagnosis Relapse}"
DEFS="${AMLGRN_DEFS:-Diagnosis Relapse}"
mkdir -p "$OUT"

partof() { local f="$1"; printf '%s.part.%s' "${f%.*}" "${f##*.}"; }

FAIL=0

# ---- [A] mudata per healthy stage, from the REPAIRED GEX -------------------------------
# This step only reads the healthy ATAC chain. The cisTopic object depends on ATAC alone, not on
# the pairing, so it is reused rather than rebuilt.
for STAGE in $STAGES; do
  MUD="$HOUT/mudata_${STAGE}_healthy_repair.h5mu"
  GEX="$HOUT/paired_rna_${STAGE}.healthy_repair.h5ad"
  CTO="$HOUT/cistopic_obj_${STAGE}.pkl"
  echo "===== [A] prepare_GEX_ACC  $STAGE ====="
  if [ -s "$MUD" ]; then echo "  already present -- skip"; continue; fi
  [ -s "$GEX" ] || { echo "  FAILED: missing $GEX"; FAIL=$((FAIL+1)); continue; }
  [ -s "$CTO" ] || { echo "  FAILED: missing $CTO -- the healthy ATAC chain is NOT present and this"
                     echo "          script refuses to rebuild an 8-day topic model implicitly";
                     FAIL=$((FAIL+1)); continue; }
  scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$CTO" \
      --GEX_anndata_fname "$GEX" \
      --out_file "$(partof "$MUD")" \
      --do_not_use_raw_for_GEX_anndata \
    || { echo "  FAILED: prepare_GEX_ACC $STAGE"; FAIL=$((FAIL+1)); continue; }
  mv -f "$(partof "$MUD")" "$MUD" || { echo "  FAILED: commit mudata $STAGE"; FAIL=$((FAIL+1)); continue; }
  echo "  -> $MUD"
done
[ "$FAIL" -eq 0 ] || { echo "==== ABORTING: mudata stage incomplete ===="; exit 1; }

# ---- [B] the cross-scored 2x2 ------------------------------------------------------------
run_one () {
  local HCELLS="$1" MDEFS="$2"
  local MUD="$HOUT/mudata_${HCELLS}_healthy_repair.h5mu"
  local EREG="$MALIG/eRegulons_${MDEFS}_repair.tsv"
  local AUC="$OUT/eRegulon_AUC_H${HCELLS}_repair_defM${MDEFS}.h5mu"

  echo "===== [B] healthy ${HCELLS} cells x malignant ${MDEFS} REPAIR defs ====="
  [ -s "$MUD" ]  || { echo "  FAILED: missing $MUD";  return 1; }
  [ -s "$EREG" ] || { echo "  FAILED: missing $EREG"; return 1; }
  if [ -s "$AUC" ]; then echo "  already present -- skip"; return 0; fi

  scenicplus grn_inference AUCell \
      --eRegulon_fname "$EREG" \
      --multiome_mudata_fname "$MUD" \
      --aucell_out_fname "$(partof "$AUC")" \
      --n_cpu 1 || { echo "  FAILED: AUCell H${HCELLS} x M${MDEFS}"; return 1; }

  mv -f "$(partof "$AUC")" "$AUC" || { echo "  FAILED: commit H${HCELLS} x M${MDEFS}"; return 1; }
  echo "  -> $AUC"
  return 0
}

# Healthy Relapse is the larger unit (34,335 repaired cells vs 25,350), so run it first so that a
# mis-sized job fails early.
for HCELLS in Relapse Diagnosis; do
  case " $STAGES " in *" $HCELLS "*) ;; *) continue ;; esac
  for MDEFS in $DEFS; do
    run_one "$HCELLS" "$MDEFS" || FAIL=$((FAIL+1))
  done
done

echo
echo "==== HEALTHY REPAIR CROSS-SCORING SUMMARY ===="
for HCELLS in $STAGES; do
  for MDEFS in $DEFS; do
    f="$OUT/eRegulon_AUC_H${HCELLS}_repair_defM${MDEFS}.h5mu"
    if [ -s "$f" ]; then printf '  OK      %s  (%s)\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
    else                 printf '  MISSING %s\n' "$(basename "$f")"; fi
  done
done
ls "$OUT"/*.part.* >/dev/null 2>&1 && echo "  NOTE: stray .part.* present (a killed run; safe to delete)"

[ "$FAIL" -eq 0 ] || { echo "==== $FAIL run(s) FAILED ===="; exit 1; }
echo "==== HEALTHY REPAIR CROSS-SCORING COMPLETE -- Gene_based is the only valid modality ===="
