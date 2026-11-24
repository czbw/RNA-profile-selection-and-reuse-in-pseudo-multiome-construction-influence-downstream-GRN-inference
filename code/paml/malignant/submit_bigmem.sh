#!/usr/bin/env bash
#SBATCH -J stage2_malig_bigmem
#SBATCH --partition=bigmem
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20    # bigmem MaxMemPerCPU=64430M -> 20 CPUs unlocks up to ~1.26 TB
#SBATCH --mem=1200G           # generous for malignant Diagnosis (~260 GB dense x 2-3x peak)
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}/paml/malignant"
SHARED="${AMLGRN_PROJECT}"            # Stage-1 top-level: holds miniforge.sif + envs + shared data
SIF="$SHARED/miniforge.sif"
mkdir -p "$ROOT/logs"

# --- pre-flight ---
[ -f "$SIF" ]                    || { echo "ERROR: $SIF not found (Stage-1 container)."; exit 1; }
[ -d "$SHARED/envs/scenicplus" ] || { echo "ERROR: $SHARED/envs/scenicplus missing (Stage-1 env)."; exit 1; }
[ -d "$SHARED/envs/macs3env" ]   || { echo "ERROR: $SHARED/envs/macs3env missing (Stage-1 env)."; exit 1; }
[ -e "$ROOT/data/fragments" ]    || { echo "ERROR: $ROOT/data not linked. Run: bash link_shared_assets.sh"; exit 1; }

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
# Keep the DAR branch's cell subsampling from firing unless the matrix will not fit.
# The budget is set just under the memory request; if Diagnosis runs out of memory, lower it
# (more subsampling) rather than raise it. Do not export this when submitting to 'long'.
export APPTAINERENV_AMLGRN_IMPUTE_BUDGET_GB=250

# Bind the SHARED parent so the bundle's symlinks (data/, databases/) resolve inside the container.
echo ">> bigmem: 20 CPUs / 1200 G, impute budget 250 GB"
echo ">> launching run_all.sh inside $SIF  (stage arg: '${1:-<all>}')"
apptainer exec --bind "$SHARED" "$SIF" bash "$ROOT/run_all.sh" "$@"
