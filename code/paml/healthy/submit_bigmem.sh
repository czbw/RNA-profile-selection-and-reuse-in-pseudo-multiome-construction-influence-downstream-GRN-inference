#!/usr/bin/env bash
#SBATCH -J armA_bigmem
#SBATCH --partition=bigmem
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16   # bigmem MaxMemPerCPU=64430M -> 16 CPUs allows up to ~1006 G
#SBATCH --mem=900G           # ~3x the projected Remission peak; keeps all 65,187 cells
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
SIF="$ROOT/miniforge.sif"

# --- pre-flight: fail loudly with a clear message rather than deep inside the run ---
[ -f "$SIF" ]                  || { echo "ERROR: $SIF not found. Run the one-time 'apptainer pull' (or scp it from WSL)."; exit 1; }
[ -d "$ROOT/envs/scenicplus" ] || { echo "ERROR: $ROOT/envs/scenicplus missing. Build the envs (see submit_container.sh)."; exit 1; }
[ -d "$ROOT/envs/macs3env" ]   || { echo "ERROR: $ROOT/envs/macs3env missing. Build the envs (see submit_container.sh)."; exit 1; }

# APPTAINERENV_* injects a variable into the container.
export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"

# Keep every cell in the DAR branch. 200 GB > Remission's 134 GB dense matrix, so region_sets.py's
# subsampling guard never fires. (Do not export this when submitting to 'long': 240 G cannot hold it.)
export APPTAINERENV_AMLGRN_IMPUTE_BUDGET_GB=200

echo ">> bigmem: 16 CPUs / 900 G, impute budget 200 GB (no cell subsampling)"
echo ">> launching run_all.sh inside $SIF  (stage arg: '${1:-<all>}')"
apptainer exec --bind "$ROOT" "$SIF" bash "$ROOT/paml/healthy/run_all.sh" "$@"
