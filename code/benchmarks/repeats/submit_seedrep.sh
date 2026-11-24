#!/bin/bash
#SBATCH --job-name=pbmc_seedrep
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
SRC=$ROOT/benchmarks/pbmc          # the completed ladder — READ ONLY
BUNDLE=$ROOT/benchmarks/repeats       # this bundle — everything new is written here
SIF=$ROOT/miniforge.sif

export TMPDIR=$ROOT/sbtmp
mkdir -p "$TMPDIR" "$BUNDLE/logs" "$BUNDLE/out"

export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"        # else EnvironmentNameNotFound
export APPTAINERENV_AMLGRN_SRC_BUNDLE="$SRC"
# chain_guard.sh resolves the chain as "$AMLGRN_ROOT/out", so AMLGRN_ROOT must point at the source
# bundle (benchmarks/pbmc); this bundle has no chain of its own.
export APPTAINERENV_AMLGRN_ROOT="$SRC"
export APPTAINERENV_AMLGRN_REP_OUT="$BUNDLE/out"
export APPTAINERENV_AMLGRN_REP_TMP="$TMPDIR"
export APPTAINERENV_AMLGRN_DB="$ROOT/data/reference/motifs"
export APPTAINERENV_AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-16}"
export APPTAINERENV_AMLGRN_STAGE="${AMLGRN_STAGE:-pbmc10k}"
export APPTAINERENV_AMLGRN_ARM="${AMLGRN_ARM:-truth}"
export APPTAINERENV_AMLGRN_REP="${AMLGRN_REP:-rep2}"
export APPTAINERENV_TMPDIR="$TMPDIR"
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$TMPDIR"

STAGE="${AMLGRN_STAGE:-pbmc10k}"
ARM="${AMLGRN_ARM:-truth}"

# Preflight on the host, before entering the container, so a missing file fails fast with a clear
# message.
for p in "$SIF" "$ROOT/envs/scenicplus" \
         "$BUNDLE/scripts/rerun_native.sh" \
         "$BUNDLE/scripts/compare_replicates.py" \
         "$SRC/scripts/chain_guard.sh"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# the four reused inputs, and the ranking database
for p in "$SRC/out/mudata_${STAGE}_${ARM}.h5mu" \
         "$SRC/out/search_space_${STAGE}_${ARM}.tsv" \
         "$SRC/out/tf_names_${STAGE}_${ARM}.txt" \
         "$SRC/out/cistromes_direct_${STAGE}_${ARM}.h5ad" \
         "$ROOT/data/reference/motifs/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"; do
  [ -e "$p" ] || { echo "MISSING REUSED INPUT: $p"; exit 2; }
done

# replicate 1 must exist, or there is nothing to compare against
for p in "$SRC/out/r2g_${STAGE}_${ARM}.tsv" "$SRC/out/tf2g_${STAGE}_${ARM}.tsv" \
         "$SRC/out/eRegulons_${STAGE}_${ARM}.tsv"; do
  [ -e "$p" ] || { echo "MISSING REPLICATE 1: $p — there is nothing to compare against"; exit 2; }
done

# chain_guard verify needs the recorded reference and has no fallback to the version-controlled
# .pin (it exits 5 without it); check on the host so the failure is named here.
CGREF="$SRC/out/atac_chain_${STAGE}.sha256"
if [ ! -s "$CGREF" ]; then
  echo "MISSING: $CGREF"
  echo "  chain_guard.sh verify requires this file and does not fall back to the .pin."
  echo "  It is written by:  AMLGRN_ROOT=$SRC bash $SRC/scripts/chain_guard.sh record $STAGE"
  echo "  Do NOT run 'record' blindly — it re-hashes whatever chain is on disk NOW. Compare the"
  echo "  result against $SRC/scripts/atac_chain_${STAGE}.sha256.pin first; the pin is the chain"
  echo "  the published arms consumed, and 'record' aborts (exit 4) if they disagree."
  exit 2
fi
echo ">> preflight OK: scripts, reused inputs, ranking db, and replicate 1 all present"

echo ">> entering $SIF once; replicate ${AMLGRN_REP:-rep2} of ${STAGE}_${ARM}"
apptainer exec --bind "$ROOT" "$SIF" bash "$BUNDLE/scripts/rerun_native.sh"
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== SEED REPLICATE COMPLETE ===="
  cat "$BUNDLE/out/seedrep_ceiling.json" 2>/dev/null
else
  echo "==== SEED REPLICATE FAILED (exit $rc) ===="
fi
exit "$rc"
