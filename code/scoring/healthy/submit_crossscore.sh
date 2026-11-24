#!/usr/bin/env bash
#SBATCH -J s3b_crossscore
#SBATCH --partition=bigmem
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    # as allocated in the earlier run (sacct AllocCPUS 16)
#SBATCH --mem=900G            # as requested in the earlier run; its sacct MaxRSS for the batch step was 99298212K
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}/scoring/healthy"
SHARED="${AMLGRN_PROJECT}"
HEALTHY="${AMLGRN_HEALTHY:-$SHARED/out}"
MALIG="$SHARED/paml/malignant/out"
SIF="$SHARED/miniforge.sif"
mkdir -p "$ROOT/logs" "$ROOT/out" "$SHARED/s3btmp"

[ -f "$SIF" ]                    || { echo "ERROR: $SIF not found."; exit 1; }
[ -d "$SHARED/envs/scenicplus" ] || { echo "ERROR: $SHARED/envs/scenicplus missing."; exit 1; }
[ -s "$ROOT/scripts/crossscore.sh" ] || { echo "ERROR: crossscore.sh not uploaded."; exit 1; }

# Refuse to run unless the pre-flight check has passed (avoids spending bigmem on an
# unverified gene universe).
V="$ROOT/out/s3b_preflight_verdict.txt"
if [ ! -s "$V" ] || ! grep -q '^PASS' "$V"; then
  echo "ERROR: pre-flight has not passed. Run 'sbatch submit_preflight.sh' first and read"
  echo "       $V . Override only with a documented reason."
  [ -s "$V" ] && { echo "       current verdict:"; sed 's/^/         /' "$V"; }
  exit 4
fi

miss=0
for f in "$HEALTHY/mudata_Diagnosis_link.h5mu" "$HEALTHY/mudata_Relapse_link.h5mu" \
         "$MALIG/eRegulons_Diagnosis_link.tsv" "$MALIG/eRegulons_Relapse_link.tsv"; do
  [ -s "$f" ] || { echo "MISSING $f"; miss=1; }
done
[ "$miss" -eq 0 ] || { echo "ERROR: inputs missing."; exit 2; }

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
export APPTAINERENV_AMLGRN_HEALTHY="$HEALTHY"
export APPTAINERENV_AMLGRN_MALIG="$MALIG"
export APPTAINERENV_AMLGRN_S3B_OUT="$ROOT/out"
export APPTAINERENV_TMPDIR="$SHARED/s3btmp"          # short path: AF_UNIX sockets cap at 107 bytes
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$SHARED/s3btmp"
export APPTAINERENV_MPLCONFIGDIR="$ROOT/out/.mpl"

echo ">> bigmem: 16 CPUs / 900 G — Stage-3B healthy-vs-malignant AUCell (4 runs)"
apptainer exec --bind "$SHARED" "$SIF" bash -c '
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate scenicplus
  bash '"$ROOT"'/scripts/crossscore.sh'
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== STAGE 3B CROSS-SCORING COMPLETE ===="
  ls -lh "$ROOT"/out/eRegulon_AUC_H*_link_defM*.h5mu 2>/dev/null
  echo "Download all four plus s3b_preflight_report.csv, then run the bystander comparison"
  echo "locally against the Stage-2/Stage-3 malignant files. Gene_based ONLY."
else
  echo "==== STAGE 3B CROSS-SCORING FAILED (exit $rc) ===="
fi
exit "$rc"
