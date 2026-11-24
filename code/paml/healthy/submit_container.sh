#!/usr/bin/env bash
#SBATCH -J armA_full
#SBATCH --partition=long
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G          # raised from 128G after an out-of-memory failure; 'long' nodes allow it
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
[ -f "$SIF" ]            || { echo "ERROR: $SIF not found. Run the one-time 'apptainer pull' (or scp it from WSL)."; exit 1; }
[ -d "$ROOT/envs/scenicplus" ] || { echo "ERROR: $ROOT/envs/scenicplus missing. Build the envs (one-time step 2)."; exit 1; }
[ -d "$ROOT/envs/macs3env" ]   || { echo "ERROR: $ROOT/envs/macs3env missing. Build the envs (one-time step 2)."; exit 1; }

# CONDA_ENVS_DIRS is injected into the container so `conda activate scenicplus` /
# `conda run -n macs3env` (inside run_all.sh) resolve to the host-built envs.
export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"

echo ">> launching run_all.sh inside $SIF  (stage arg: '${1:-<all>}')"
apptainer exec --bind "$ROOT" "$SIF" bash "$ROOT/paml/healthy/run_all.sh" "$@"
