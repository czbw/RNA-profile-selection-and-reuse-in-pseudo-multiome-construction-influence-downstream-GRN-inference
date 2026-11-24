#!/bin/bash
#SBATCH --job-name=aml_repair
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=20
#SBATCH --mem=1200G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
STAGE2=$ROOT/paml/malignant      # owns the shared ATAC chain + out/
BUNDLE=$ROOT/paml/corrected            # owns the repair scripts + the repaired cell_data
SIF=$ROOT/miniforge.sif

export TMPDIR=$ROOT/pbtmp
mkdir -p "$TMPDIR" "$BUNDLE/logs"

export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"          # else EnvironmentNameNotFound
export APPTAINERENV_AMLGRN_ROOT="$STAGE2"                 # chain + out/ live here
export APPTAINERENV_AMLGRN_REPAIR_BUNDLE="$BUNDLE"
export APPTAINERENV_AMLGRN_CONFIG="config.yaml"
export APPTAINERENV_AMLGRN_DB="$ROOT/data/reference/motifs"
export APPTAINERENV_AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-20}"
export APPTAINERENV_AMLGRN_STAGES="${*:-Diagnosis Relapse}"
export APPTAINERENV_TMPDIR="$TMPDIR"
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$TMPDIR"

# Preflight: the builder imports `common` and reads config.yaml from its own directory, so
# check that these files exist before entering the container and fail fast with a clear message.
for p in "$SIF" "$ROOT/envs/scenicplus" "$BUNDLE/run_corrected.sh" \
         "$BUNDLE/scripts/run_grn.sh" "$BUNDLE/scripts/build_paired_rna.py" \
         "$BUNDLE/scripts/common.py" "$BUNDLE/scripts/config.yaml"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# the chain itself, and the uploaded assignment -- both are hard requirements
for STAGE in ${*:-Diagnosis Relapse}; do
  for p in "$STAGE2/out/cistopic_obj_${STAGE}.pkl" "$STAGE2/out/menr_cistarget_${STAGE}.hdf5" \
           "$STAGE2/out/region_sets_${STAGE}" \
           "$BUNDLE/uploads/cell_data_${STAGE}.repaired.csv"; do
    [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
  done
done
echo ">> preflight OK: scripts, dependencies, ATAC chain and assignment tables all present"

echo ">> entering $SIF once; stages: ${*:-Diagnosis Relapse}"
apptainer exec --bind "$ROOT" "$SIF" bash "$BUNDLE/run_corrected.sh"
