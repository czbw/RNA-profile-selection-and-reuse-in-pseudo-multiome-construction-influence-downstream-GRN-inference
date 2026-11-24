#!/usr/bin/env bash
#SBATCH -J stage2_malig_long
#SBATCH --partition=long
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G           # 'long' ceiling (~288 GB physical). Sized for Relapse.
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}/paml/malignant"
SHARED="${AMLGRN_PROJECT}"
SIF="$SHARED/miniforge.sif"
mkdir -p "$ROOT/logs"

# --- pre-flight: partition/stage sizing ---
# No argument means run_all.sh runs both stages, including Diagnosis; refuse that as well.
STAGE_ARG="${1:-}"
if [ -z "$STAGE_ARG" ] || [ "$STAGE_ARG" = "Diagnosis" ]; then
  if [ "${AMLGRN_ALLOW_LONG_DIAGNOSIS:-0}" != "1" ]; then
    echo "ERROR: malignant Diagnosis does not fit on 'long' (240 G) — a previous run was OOM-killed"
    echo "       in AUCell here. Submit it to bigmem instead:"
    echo "         sbatch submit_bigmem.sh Diagnosis"
    echo "       (Relapse is what this script is for: sbatch submit_container.sh Relapse)"
    echo "       To override deliberately: --export=ALL,AMLGRN_ALLOW_LONG_DIAGNOSIS=1"
    exit 4
  fi
  echo ">> WARNING: AMLGRN_ALLOW_LONG_DIAGNOSIS=1 — running Diagnosis on 'long' against advice."
fi

[ -f "$SIF" ]                    || { echo "ERROR: $SIF not found (Stage-1 container)."; exit 1; }
[ -d "$SHARED/envs/scenicplus" ] || { echo "ERROR: $SHARED/envs/scenicplus missing (Stage-1 env)."; exit 1; }
[ -d "$SHARED/envs/macs3env" ]   || { echo "ERROR: $SHARED/envs/macs3env missing (Stage-1 env)."; exit 1; }
[ -e "$ROOT/data/fragments" ]    || { echo "ERROR: $ROOT/data not linked. Run: bash link_shared_assets.sh"; exit 1; }

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
# Do not export AMLGRN_IMPUTE_BUDGET_GB here: the config default (64 GB) lets the DAR branch
# subsample. This protects region_sets.py only; see the AUCell note in scripts/run_grn.sh.

echo ">> long: 16 CPUs / 240 G (Relapse-sized; Diagnosis is refused above)"
echo ">> launching run_all.sh inside $SIF  (stage arg: '${1:-<all>}')"
apptainer exec --bind "$SHARED" "$SIF" bash "$ROOT/run_all.sh" "$@"
