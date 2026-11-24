#!/usr/bin/env bash
#SBATCH -J s3_crossscore
#SBATCH --partition=bigmem
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20    # bigmem MaxMemPerCPU=64430M -> 20 CPUs unlocks up to ~1.26 TB
#SBATCH --mem=1200G
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}/scoring/published"
SHARED="${AMLGRN_PROJECT}"
S2="$SHARED/paml/malignant/out"
SIF="$SHARED/miniforge.sif"
mkdir -p "$ROOT/logs" "$ROOT/out" "$SHARED/s3tmp"

[ -f "$SIF" ]                    || { echo "ERROR: $SIF not found."; exit 1; }
[ -d "$SHARED/envs/scenicplus" ] || { echo "ERROR: $SHARED/envs/scenicplus missing."; exit 1; }
[ -s "$ROOT/scripts/crossscore.sh" ] || { echo "ERROR: crossscore.sh not uploaded."; exit 1; }

# Refuse to run unless the pre-flight check has passed (avoids spending 1.2 TB of memory on an
# unverified gene universe).
V="$ROOT/out/s3_preflight_verdict.txt"
if [ ! -s "$V" ] || ! grep -q '^PASS' "$V"; then
  echo "ERROR: pre-flight has not passed. Run 'sbatch submit_preflight.sh' first and read"
  echo "       $V . Override only with a documented reason."
  [ -s "$V" ] && { echo "       current verdict:"; sed 's/^/         /' "$V"; }
  exit 4
fi

miss=0
for f in "$S2/mudata_Diagnosis_link.h5mu" "$S2/mudata_Relapse_link.h5mu" \
         "$S2/eRegulons_Diagnosis_link.tsv" "$S2/eRegulons_Relapse_link.tsv"; do
  [ -s "$f" ] || { echo "MISSING $f"; miss=1; }
done
[ "$miss" -eq 0 ] || { echo "ERROR: Stage-2 inputs missing."; exit 2; }

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
export APPTAINERENV_AMLGRN_S2="$S2"
export APPTAINERENV_AMLGRN_S3_OUT="$ROOT/out"
export APPTAINERENV_TMPDIR="$SHARED/s3tmp"          # short path: AF_UNIX sockets cap at 107 bytes
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$SHARED/s3tmp"
export APPTAINERENV_MPLCONFIGDIR="$ROOT/out/.mpl"

echo ">> bigmem: 20 CPUs / 1200 G — Stage-3 cross-stage AUCell (2 runs)"
apptainer exec --bind "$SHARED" "$SIF" bash -c '
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate scenicplus
  bash '"$ROOT"'/scripts/crossscore.sh'
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== STAGE 3 CROSS-SCORING COMPLETE ===="
  ls -lh "$ROOT"/out/eRegulon_AUC_*_link_def*.h5mu 2>/dev/null
  echo "Download these two, then fit the Stage-3 model locally. Gene_based ONLY."
else
  echo "==== STAGE 3 CROSS-SCORING FAILED (exit $rc) ===="
fi
exit "$rc"
