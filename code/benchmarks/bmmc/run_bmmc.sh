#!/usr/bin/env bash
# Driver for the BMMC ladder: build the shared chromatin chain once, freeze it, then run one GRN
# per arm against the frozen chain. The ladder's claim is that any difference between arms is due
# to the cell pairing alone, which requires every arm to consume the identical chain.
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus

BUNDLE="${AMLGRN_ROOT:?submit_bmmc.sh must export AMLGRN_ROOT}"
STAGE="${AMLGRN_STAGE:-bmmc}"
CFG="${AMLGRN_CONFIG:-config.yaml}"
SEED="${AMLGRN_SEED:-555}"
ARMS="${AMLGRN_ARMS:-truth rand_celltype_1to1 rand_celltype_donor_1to1 rand_global}"
NCPU="${AMLGRN_NCPU:-20}"
DB="${AMLGRN_DB:-$BUNDLE/databases}"
OUT="$BUNDLE/out"; WORK="$BUNDLE/work"; UP="$BUNDLE/uploads"; DATA="$BUNDLE/data"
TMP="$WORK/grn_tmp_${STAGE}_shared"
# optional early exit, so the chain can be built in its own allocation and the arms submitted
# afterwards against a frozen chain:  AMLGRN_STOP_AFTER=chain sbatch submit_bmmc.sh
STOP_AFTER="${AMLGRN_STOP_AFTER:-}"

RANK="$DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
MOTIF_TBL="$DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

mkdir -p "$OUT" "$WORK" "$TMP"
export AMLGRN_ROOT="$BUNDLE" AMLGRN_STAGE="$STAGE" AMLGRN_CONFIG="$CFG" \
       AMLGRN_NCPU="$NCPU" AMLGRN_DB="$DB" AMLGRN_SEED="$SEED"

echo "=================================================================="
echo " BMMC LADDER   $(date)"
echo " python  : $(which python)"
echo " bundle  : $BUNDLE"
echo " stage   : $STAGE   seed: $SEED   ncpu: $NCPU"
echo " arms    : $ARMS"
echo " DB      : $DB"
echo " TMPDIR  : ${TMPDIR:-unset}"
echo " dataset : GSE194122 BMMC 10x Multiome — 69,249 cells, 13,431 GEX features,"
echo "           116,490 ATAC peaks, 13 batches, 10 donors, 4 sites, 22 cell types"
echo "=================================================================="

FAIL=0
IMPOSSIBLE=""

# Calling convention (required): cd into scripts/ and invoke bare filenames. The python steps resolve
# sibling modules and the config relative to their own directory, and run_grn_bmmc.sh invokes helpers
# without a directory prefix; any other cwd breaks an arm hours in.
cd "$BUNDLE/scripts" || { echo "FATAL: no $BUNDLE/scripts"; exit 3; }

# `truth` is the barcode-identity pairing, available only from a real multiome, so it has no uploaded
# assignment table. Metacell arms have no per-cell assignment: SCENIC+ builds the metacells itself
# inside prepare_GEX_ACC.
needs_assignment() {
  case "$1" in
    truth) return 1 ;;
    pseudomultiome|pseudomultiome_k*) return 1 ;;
    *) return 0 ;;
  esac
}
is_metacell() {
  case "$1" in pseudomultiome|pseudomultiome_k*) return 0 ;; *) return 1 ;; esac
}
assignment_csv() { echo "$UP/arm_${1}_${STAGE}_seed${SEED}_full.csv"; }
arm_rna()        { echo "$OUT/paired_rna_${STAGE}.${1}.h5ad"; }

partof() { local f="$1"; printf '%s.part.%s' "${f%.*}" "${f##*.}"; }

# [1/4] ASSERT. Nothing is written in this phase.
echo; echo "---- [1/4] asserting inputs (writes nothing) ----"

# the scripts this run will actually execute
for f in "$BUNDLE/scripts/$CFG" \
         "$BUNDLE/scripts/atac_chain.py" \
         "$BUNDLE/scripts/region_sets.py" \
         "$BUNDLE/scripts/build_arm_rna.py" \
         "$BUNDLE/scripts/run_grn_bmmc.sh" \
         "$BUNDLE/scripts/chain_guard.sh"; do
  [ -s "$f" ] || { echo "  MISSING SCRIPT $f"; FAIL=$((FAIL+1)); }
done

# the source matrix and the join table
# b00_cell_metadata.csv carries the verified barcode remap for all 69,249 cells:
#   h5ad AAACAGCCAATTAAGG-1-s1d1  ->  fragment AAACAGCCAATTAAGG-1
# the stripped suffix is the batch. Every downstream join keys on that column, not on obs/_index.
H5AD="$DATA/multiome_BMMC.h5ad"
[ -s "$H5AD" ] || { echo "  MISSING $H5AD (gzip -dc the GSE194122 .h5ad.gz on CUBIC; ~2.90 GiB)"; FAIL=$((FAIL+1)); }
[ -s "$UP/b00_cell_metadata.csv" ] || { echo "  MISSING $UP/b00_cell_metadata.csv (scp it from benchmarks/bmmc/out/)"; FAIL=$((FAIL+1)); }

# the cisTarget databases. BMMC is GRCh38, so the feathers already on CUBIC are reused as-is;
# nothing genome-specific is staged for this bundle.
for f in "$RANK" "$MOTIF_TBL" "$DB/hg38_genome_annotation.tsv" "$DB/hg38_chromsizes_3col.tsv"; do
  [ -s "$f" ] || { echo "  MISSING DATABASE $f"; FAIL=$((FAIL+1)); }
done

# every requested arm, checked NOW rather than when the loop reaches it
for ARM in $ARMS; do
  if needs_assignment "$ARM"; then
    CSV="$(assignment_csv "$ARM")"
    if [ -s "$CSV" ]; then
      echo "  [$ARM] assignment: $(basename "$CSV") ($(wc -l < "$CSV") lines)"
    else
      echo "  [$ARM] MISSING $CSV -- build it locally and scp it to uploads/"
      FAIL=$((FAIL+1))
    fi
  elif is_metacell "$ARM"; then
    G="$OUT/native_rna_${STAGE}.metacell.h5ad"
    [ -s "$G" ] && echo "  [$ARM] metacell GEX present (no assignment table by design)" \
                || { echo "  [$ARM] MISSING $G"; FAIL=$((FAIL+1)); }
  else
    echo "  [$ARM] barcode-identity pairing -- no assignment table by design"
  fi
done

[ "$FAIL" -eq 0 ] || { echo; echo "==== ABORTING: $FAIL input(s) missing -- nothing was written ===="; exit 1; }
echo "  all inputs present for: $ARMS"

# [2/4] Shared chromatin chain: built once, hashed, never rebuilt.
# It is ATAC-only, hence pairing-independent and identical for every arm. Resumable: each sub-step
# is skipped if its output already exists non-empty.
echo; echo "---- [2/4] shared chromatin chain (ATAC-only; built once for ALL arms) ----"

CIS="$OUT/cistopic_obj_${STAGE}.pkl"
RGS="$OUT/region_sets_${STAGE}"
MENR="$OUT/menr_cistarget_${STAGE}.hdf5"

# [2a] peak x cell count matrix -> cisTopic object -> LDA topic models.
# atac_chain.py reads layers/counts (not X, which is library-size-normalised + log1p and must
# not be summed) and rewrites peak names from chr1-9776-10668 to chr1:9776-10668, the only form
# pycisTopic parses.
if [ -s "$CIS" ]; then
  echo "  [2a] cistopic_obj present -- skip"
else
  echo "  [2a] atac_chain.py (peak matrix -> cisTopic -> topics)"
  t0=$(date +%s)
  python atac_chain.py --h5ad "$H5AD" --metadata "$UP/b00_cell_metadata.csv" \
    || { echo "  b_01_atac_chain FAILED"; echo "==== ABORTING ===="; exit 1; }
  [ -s "$CIS" ] || { echo "  b_01 returned 0 but $CIS is missing/empty"; echo "==== ABORTING ===="; exit 1; }
  echo "  [2a] took $(( ($(date +%s) - t0) / 60 )) min"
fi

# [2b] binarised topics + per-cell-type DARs
if [ -d "$RGS" ]; then
  echo "  [2b] region_sets present -- skip"
else
  echo "  [2b] region_sets.py"
  t0=$(date +%s)
  python region_sets.py \
    || { echo "  b_02_region_sets FAILED"; echo "==== ABORTING ===="; exit 1; }
  [ -d "$RGS" ] || { echo "  b_02 returned 0 but $RGS is missing"; echo "==== ABORTING ===="; exit 1; }
  echo "  [2b] took $(( ($(date +%s) - t0) / 60 )) min"
fi

# [2c] motif enrichment, hoisted out of the arm loop. Written to a .part sidecar and moved into place
# only on exit 0, so a walltime/OOM kill mid-write cannot leave a truncated .hdf5 that a later
# resubmit would accept as complete and hash into the frozen chain.
if [ -s "$MENR" ]; then
  echo "  [2c] motif enrichment present -- skip"
else
  echo "  [2c] motif_enrichment_cistarget (ONCE, shared by every arm)"
  t0=$(date +%s)
  scenicplus grn_inference motif_enrichment_cistarget \
    --region_set_folder "$RGS" \
    --cistarget_db_fname "$RANK" \
    --output_fname_cistarget_result "$(partof "$MENR")" \
    --temp_dir "$TMP" --species homo_sapiens \
    --path_to_motif_annotations "$MOTIF_TBL" \
    --annotation_version v10nr_clust --n_cpu 1 \
    || { echo "  motif_enrichment_cistarget FAILED"; rm -f "$(partof "$MENR")"; echo "==== ABORTING ===="; exit 1; }
  mv -f "$(partof "$MENR")" "$MENR" \
    || { echo "  motif_enrichment_cistarget (commit) FAILED"; echo "==== ABORTING ===="; exit 1; }
  echo "  [2c] took $(( ($(date +%s) - t0) / 60 )) min"
fi

# [2d] Freeze the chain. Hash recorded before the first arm; run_grn_bmmc.sh re-verifies at every
# arm entry and this script re-verifies after the last arm.
if [ -s "$OUT/atac_chain_${STAGE}.sha256" ]; then
  echo "  [2d] chain hashes already recorded -- verifying instead"
  bash chain_guard.sh verify "$STAGE" "pre-arms" \
    || { echo "==== ABORTING: the recorded chain no longer matches what is on disk ===="; exit 1; }
else
  echo "  [2d] recording ATAC-chain hashes"
  bash chain_guard.sh record "$STAGE" \
    || { echo "  chain_guard record FAILED"; echo "==== ABORTING ===="; exit 1; }
fi

if [ "$STOP_AFTER" = "chain" ]; then
  echo; echo "==== AMLGRN_STOP_AFTER=chain -- chain built and frozen, stopping before the arms ===="
  exit 0
fi

# [3/4] PER-ARM RNA OBJECTS — all of them, before any GRN runs.
# BMMC is a true multiome, so a rebuilt pairing is a permutation within a stratum: n_ATAC == n_RNA,
# a bijection is always feasible, and there is no scarcity floor (unlike an unpaired cohort where
# a stratum can run out of donors). If b_03 reports dropped cells, the assignment table is wrong.
echo; echo "---- [3/4] per-arm RNA objects ----"
BUILT=0
for ARM in $ARMS; do
  if is_metacell "$ARM"; then
    echo "  [$ARM] metacell arm -- SCENIC+ builds the metacells inside prepare_GEX_ACC; nothing to do"
    continue
  fi
  RNA="$(arm_rna "$ARM")"
  if [ -s "$RNA" ]; then
    echo "  [$ARM] $(basename "$RNA") present -- skip"
    continue
  fi
  if needs_assignment "$ARM"; then
    python build_arm_rna.py --stage "$STAGE" --arm "$ARM" --seed "$SEED" \
      --h5ad "$H5AD" --metadata "$UP/b00_cell_metadata.csv" \
      --assignment "$(assignment_csv "$ARM")" \
      || { echo "  [$ARM] b_03_build_arm_rna FAILED"; FAIL=$((FAIL+1)); continue; }
  else
    python build_arm_rna.py --stage "$STAGE" --arm "$ARM" --seed "$SEED" \
      --h5ad "$H5AD" --metadata "$UP/b00_cell_metadata.csv" \
      || { echo "  [$ARM] b_03_build_arm_rna FAILED"; FAIL=$((FAIL+1)); continue; }
  fi
  [ -s "$RNA" ] || { echo "  [$ARM] b_03 returned 0 but $(basename "$RNA") is missing/empty"; FAIL=$((FAIL+1)); continue; }
  BUILT=$((BUILT+1))
done
[ "$FAIL" -eq 0 ] || { echo; echo "==== ABORTING: $FAIL arm(s) could not be built -- no GRN was started ===="; exit 1; }
echo "  built $BUILT new arm object(s); every requested arm has an RNA object"

# [4/4] GRN HALF, PER ARM — independent and resumable.
# run_grn_bmmc.sh re-verifies the frozen chain at entry and aborts non-zero on any mismatch. An arm is
# skipped if its AUCell output already exists non-empty, so a resubmit resumes at the first unfinished arm.
echo; echo "---- [4/4] GRN half, per arm ----"
for ARM in $ARMS; do
  AUC="$OUT/eRegulon_AUC_${STAGE}_${ARM}.h5mu"
  if [ -s "$AUC" ]; then
    echo "  [$ARM] already complete -- skip"
    continue
  fi
  t0=$(date +%s)
  bash run_grn_bmmc.sh "$STAGE" "$ARM"
  rc=$?
  if [ "$rc" -eq 45 ]; then
    echo "  [$ARM] eGRN structurally impossible -- recorded, not counted as a failure"
    IMPOSSIBLE="$IMPOSSIBLE $ARM"
  elif [ "$rc" -ne 0 ]; then
    echo "  [$ARM] GRN FAILED (rc=$rc)"; FAIL=$((FAIL+1)); continue
  fi
  echo "  [$ARM] GRN half took $(( ($(date +%s) - t0) / 60 )) min"
done

# ---- re-verify the chain AFTER the last arm ---------------------------------------------------
# This is the other half of the invariant: recording the chain hash before the arms shows what they
# started from; verifying it after shows nothing changed underneath them.
echo; echo "---- re-verifying the chromatin chain after the last arm ----"
bash chain_guard.sh verify "$STAGE" "post-arms" || FAIL=$((FAIL+1))

echo
echo "==== BMMC LADDER SUMMARY  $(date) ===="
for ARM in $ARMS; do
  f="$OUT/eRegulon_AUC_${STAGE}_${ARM}.h5mu"
  if [ -s "$f" ]; then
    printf '  OK      %-40s %s\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else
    case " $IMPOSSIBLE " in
      *" $ARM "*) printf '  eGRN-IMPOSSIBLE  %s (recorded result, not a failure)\n' "$ARM" ;;
      *)          printf '  MISSING %s\n' "$(basename "$f")"; FAIL=$((FAIL+1)) ;;
    esac
  fi
done
if [ "$FAIL" -ne 0 ]; then echo "  $FAIL step(s) FAILED"; echo "==== INCOMPLETE ===="; exit 1; fi
if [ -n "$IMPOSSIBLE" ]; then
  echo "  eGRN structurally impossible (recorded as a result):$IMPOSSIBLE"
fi
echo "==== BMMC LADDER COMPLETE ===="
