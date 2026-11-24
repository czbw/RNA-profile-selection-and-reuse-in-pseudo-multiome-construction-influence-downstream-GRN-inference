#!/usr/bin/env bash
#SBATCH -J armA_setup
#SBATCH --partition=long
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
SIF="$ROOT/miniforge.sif"

# keep apptainer's cache + build tmp off the (small) default locations
export APPTAINER_CACHEDIR="$ROOT/.apptainer_cache"
export APPTAINER_TMPDIR="$ROOT/.apptainer_tmp"
# make conda-in-container put/find named envs on host disk, and write the package cache to a
# writable host dir (the container's /opt/conda/pkgs is read-only inside the SIF).
export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"
export APPTAINERENV_CONDA_PKGS_DIRS="$ROOT/.conda_pkgs"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$ROOT/envs" "$ROOT/.conda_pkgs" "$ROOT/logs"

# --- 1) pull the modern-glibc base image (idempotent) ---
if [ ! -f "$SIF" ]; then
  echo ">>>> [1] pulling miniforge base image (mksquashfs runs here, not on login)"
  apptainer pull "$SIF" docker://condaforge/miniforge3:latest
else
  echo ">>>> [1] $SIF already present — skipping pull"
fi

# --- 2) scenicplus env (inside the container, onto host envs/) ---
# Abort on failure rather than continuing to macs3/verify with a missing env
# (a solve error goes to the .err log and would otherwise go unnoticed).
if [ ! -d "$ROOT/envs/scenicplus" ]; then
  echo ">>>> [2] building 'scenicplus' env inside the container"
  apptainer exec --bind "$ROOT" "$SIF" mamba env create -n scenicplus -f "$ROOT/paml/healthy/env/environment.yml" \
    || { echo "ERROR: 'scenicplus' env build FAILED — see logs/armA_setup_*.err. Aborting."; exit 1; }
  [ -d "$ROOT/envs/scenicplus" ] \
    || { echo "ERROR: 'scenicplus' prefix not created despite no error. Aborting."; exit 1; }
else
  echo ">>>> [2] envs/scenicplus exists — skipping"
fi

if [ ! -d "$ROOT/envs/macs3env" ]; then
  echo ">>>> [3] building 'macs3env' env inside the container"
  apptainer exec --bind "$ROOT" "$SIF" mamba env create -n macs3env -f "$ROOT/paml/healthy/env/environment.macs3.yml" \
    || { echo "ERROR: 'macs3env' env build FAILED — see logs/armA_setup_*.err. Aborting."; exit 1; }
  [ -d "$ROOT/envs/macs3env" ] \
    || { echo "ERROR: 'macs3env' prefix not created despite no error. Aborting."; exit 1; }
else
  echo ">>>> [3] envs/macs3env exists — skipping"
fi

# --- 4) verify (import the core stack + locate macs3) ---
echo ">>>> [4] verifying scenicplus stack"
apptainer exec --bind "$ROOT" "$SIF" conda run -n scenicplus \
  python -c "import scenicplus, pycisTopic, pyscenic, mudata, polars; print('scenicplus', scenicplus.__version__)"
echo ">> macs3 at: $(apptainer exec --bind "$ROOT" "$SIF" conda run -n macs3env which macs3)"
echo "==== CONTAINER SETUP COMPLETE — now: sbatch submit_container.sh ===="
