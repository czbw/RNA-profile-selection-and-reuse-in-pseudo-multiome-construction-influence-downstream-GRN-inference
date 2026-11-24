#!/usr/bin/env bash
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --job-name=stage2_malig
#SBATCH --output=%x_%j.log
#SBATCH --partition=__TODO_PARTITION__
#SBATCH --account=__TODO_ACCOUNT__
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

export AMLGRN_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export AMLGRN_CONFIG=config.yaml
export AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-16}"

# Keep joblib/loky memmaps and temp files on disk, not the container's small /dev/shm.
# TMPDIR must also be short: pycisTopic's peak_calling starts ray, which creates an AF_UNIX
# socket at $TMPDIR/ray/session_<ts>_<pid>/sockets/plasma_store, and AF_UNIX paths cannot exceed
# 107 bytes. The deep bundle path $AMLGRN_ROOT/work/tmp overflows this limit, so TMPDIR is placed
# at the shared parent directory (on disk and inside the container's bound directory).
export TMPDIR="${AMLGRN_ROOT%/*}/s2tmp"          # = /path/to/AML-GRN/s2tmp (short, bound, on disk)
export JOBLIB_TEMP_FOLDER="$TMPDIR"
export MPLCONFIGDIR="$AMLGRN_ROOT/work/.mpl"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$AMLGRN_ROOT/work"

cd "$AMLGRN_ROOT/scripts"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus
export AMLGRN_MACS="$(conda run -n macs3env which macs3)"
echo ">> AMLGRN_ROOT=$AMLGRN_ROOT  ncpu=$AMLGRN_NCPU  macs3=$AMLGRN_MACS  (STAGE 2 / malignant)"

ONLY_STAGE="${1:-}"
if [ -n "$ONLY_STAGE" ]; then STAGES=("$ONLY_STAGE"); else STAGES=(Diagnosis Relapse); fi
PAIRINGS=(link pseudomultiome)

# --- preflight: malignant RNA h5ads must be staged (built by build_malignant_h5ad.py) ---
missing_rna=0
for s in Diagnosis Relapse; do
  f="$AMLGRN_ROOT/data/rna/malignant_$(echo "$s" | tr '[:upper:]' '[:lower:]').h5ad"
  [ -f "$f" ] || { echo "ERROR: missing $f — run build_malignant_h5ad.py locally and scp it (see README)"; missing_rna=1; }
done
[ "$missing_rna" -eq 0 ] || exit 2

# --- [0] fragments. Stage 1 already converted every sample into the SHARED data/fragments
#     (symlinked in). Skip re-extraction of the 50 GB tar if the converted+indexed files exist;
#     only fall back to prepare_fragments.sh if the shared fragments are somehow absent. ---
if ls "$AMLGRN_ROOT"/data/fragments/*.fragments.tsv.gz.tbi >/dev/null 2>&1; then
  echo ">>>> [0] fragments: reusing shared data/fragments ($(ls "$AMLGRN_ROOT"/data/fragments/*.tbi | wc -l) samples indexed) — skip"
else
  echo ">>>> [0] fragments: shared set missing -> extracting/converting"; bash prepare_fragments.sh
fi
echo ">>>> [1] cell_data per stage (MALIGNANT filter)"
if [ -s "$AMLGRN_ROOT/work/cell_data_Diagnosis.csv" ] && [ -s "$AMLGRN_ROOT/work/cell_data_Relapse.csv" ]; then
  echo "[stage2] cell_data present for both stages — skip"
else
  python make_cell_data.py || { echo "FATAL: make_cell_data.py failed — aborting"; exit 3; }
fi

# GEX builds fail fast: both pairings need them, so a failure here aborts the whole run
# (rather than running the ATAC chain and failing later in run_grn.sh on the missing GEX).
echo ">>>> [2] GEX builds (A paired + B native)"
gex_ok=1
for s in Diagnosis Relapse; do
  [ -s "$AMLGRN_ROOT/out/paired_rna_$s.h5ad" ] && [ -s "$AMLGRN_ROOT/out/native_rna_$s.h5ad" ] || gex_ok=0
done
if [ "$gex_ok" -eq 1 ]; then
  echo "[stage2] GEX h5ads present for both stages (paired + native) — skip"
else
  python build_paired_rna.py || { echo "FATAL: build_paired_rna.py failed (Option A GEX) — aborting"; exit 3; }
  python build_native_rna.py || { echo "FATAL: build_native_rna.py failed (Option B GEX) — aborting"; exit 3; }
fi

FAILURES=()
for S in "${STAGES[@]}"; do
  export AMLGRN_STAGE="$S"
  echo ">>>> [3/$S] ATAC chain (pseudobulk->macs3->consensus->cisTopic->topics)"
  if ! python atac_chain.py; then
    echo "SKIP [$S]: ATAC chain failed -> skipping region sets + GRN for this stage"
    FAILURES+=("$S: ATAC chain"); continue
  fi
  echo ">>>> [4/$S] region sets"
  if ! python region_sets.py; then
    echo "SKIP [$S]: region sets failed -> skipping GRN for this stage"
    FAILURES+=("$S: region sets"); continue
  fi
  for P in "${PAIRINGS[@]}"; do
    echo ">>>> [5/$S/$P] GRN half"
    if ! bash run_grn.sh "$S" "$P"; then
      echo "SKIP [$S/$P]: GRN half failed -> continuing to next pairing"
      FAILURES+=("$S/$P: GRN half")
    fi
  done
done


echo
echo "==== STAGE 2 (malignant) SUMMARY ===="
for S in "${STAGES[@]}"; do
  for P in "${PAIRINGS[@]}"; do
    F="$AMLGRN_ROOT/out/eRegulons_${S}_${P}.tsv"
    if [ -s "$F" ]; then printf '  OK      %-12s %-16s %s eRegulon rows\n' "$S" "$P" "$(($(wc -l < "$F") - 1))"
    else                 printf '  MISSING %-12s %-16s\n' "$S" "$P"; fi
  done
done
if [ ${#FAILURES[@]} -gt 0 ]; then
  echo
  echo "==== STAGE 2 INCOMPLETE — ${#FAILURES[@]} step(s) failed ===="
  for f in "${FAILURES[@]}"; do echo "  FAILED: $f"; done
  exit 1
fi
echo "==== STAGE 2 MALIGNANT COMPLETE (Dx + Rel x both pairings) ===="
