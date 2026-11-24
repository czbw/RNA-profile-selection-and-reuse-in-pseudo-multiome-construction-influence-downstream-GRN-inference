#!/usr/bin/env bash
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --job-name=armA_full
#SBATCH --output=%x_%j.log
#SBATCH --partition=__TODO_PARTITION__      # <-- CUBIC partition/QOS
#SBATCH --account=__TODO_ACCOUNT__          # <-- your allocation
#SBATCH --time=5-00:00:00                   # <-- the Remission topic model is the limiting step; size generously
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G                           # peak: Remission cisTopic + GBM region/TF_to_gene
set -uo pipefail

export AMLGRN_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export AMLGRN_CONFIG=config.yaml
export AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-16}"

# Keep joblib/loky memmaps and temp files on disk: the container's /dev/shm is small, and tmpfs
# usage counts against the SLURM memory limit. Point every temp mechanism at work/tmp.
export TMPDIR="$AMLGRN_ROOT/work/tmp"
export JOBLIB_TEMP_FOLDER="$TMPDIR"
export MPLCONFIGDIR="$AMLGRN_ROOT/work/.mpl"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR"

cd "$AMLGRN_ROOT/scripts"

# --- activate envs ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus
export AMLGRN_MACS="$(conda run -n macs3env which macs3)"
echo ">> AMLGRN_ROOT=$AMLGRN_ROOT  ncpu=$AMLGRN_NCPU  macs3=$AMLGRN_MACS"

ONLY_STAGE="${1:-}"                          # optional single-stage (array) mode
if [ -n "$ONLY_STAGE" ]; then STAGES=("$ONLY_STAGE"); else STAGES=(Diagnosis Remission Relapse); fi
PAIRINGS=(link pseudomultiome)

# --- shared, once (idempotent) ---
echo ">>>> [0] fragments (extract + convert)"; bash prepare_fragments.sh
echo ">>>> [1] cell_data per stage";           python make_cell_data.py
echo ">>>> [2] GEX builds (A paired + B native)"; python build_paired_rna.py; python build_native_rna.py

# --- per stage: ATAC chain -> region sets -> both pairings' GRN ---
# Fail fast within a stage (dependency chain ATAC -> region sets -> GRN): if an upstream step
# fails, skip its downstream steps and move to the next stage. Independent stages still run.
FAILURES=()   # every skip is recorded, summarised at the end, and the job exits non-zero

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
echo "==== ARM A SUMMARY ===="
for S in "${STAGES[@]}"; do
  for P in "${PAIRINGS[@]}"; do
    F="$AMLGRN_ROOT/out/eRegulons_${S}_${P}.tsv"
    if [ -s "$F" ]; then printf '  OK      %-12s %-16s %s eRegulon rows\n' "$S" "$P" "$(($(wc -l < "$F") - 1))"
    else                 printf '  MISSING %-12s %-16s\n' "$S" "$P"; fi
  done
done
if [ ${#FAILURES[@]} -gt 0 ]; then
  echo
  echo "==== ARM A INCOMPLETE — ${#FAILURES[@]} step(s) failed ===="
  for f in "${FAILURES[@]}"; do echo "  FAILED: $f"; done
  exit 1
fi
echo "==== ARM A FULL COMPLETE (all stages x pairings) ===="
