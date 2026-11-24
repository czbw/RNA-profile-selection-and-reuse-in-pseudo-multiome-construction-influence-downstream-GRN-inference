#!/usr/bin/env bash
#SBATCH -J s3r_preflight
#SBATCH --partition=bioinformatics
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}/scoring/corrected"
SHARED="${AMLGRN_PROJECT}"
S2="$SHARED/paml/malignant/out"
SIF="$SHARED/miniforge.sif"
ARM=repair
mkdir -p "$ROOT/logs" "$ROOT/out"

[ -f "$SIF" ]                    || { echo "ERROR: $SIF not found."; exit 1; }
[ -d "$SHARED/envs/scenicplus" ] || { echo "ERROR: $SHARED/envs/scenicplus missing."; exit 1; }
[ -s "$ROOT/scripts/preflight.py" ] || { echo "ERROR: preflight.py not uploaded."; exit 1; }

miss=0
for f in "$S2/mudata_Diagnosis_${ARM}.h5mu" "$S2/mudata_Relapse_${ARM}.h5mu" \
         "$S2/eRegulons_Diagnosis_${ARM}.tsv" "$S2/eRegulons_Relapse_${ARM}.tsv"; do
  if [ -s "$f" ]; then printf '  OK      %-46s %s\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else printf '  MISSING %s\n' "$f"; miss=1; fi
done
[ "$miss" -eq 0 ] || { echo "ERROR: Stage-2 repair inputs missing."; exit 2; }

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
export APPTAINERENV_AMLGRN_S2="$S2"
export APPTAINERENV_AMLGRN_S3_OUT="$ROOT/out"
export APPTAINERENV_AMLGRN_S3_PAIRING="$ARM"
export APPTAINERENV_MPLCONFIGDIR="$ROOT/out/.mpl"

echo ">> bioinformatics: 2 CPUs / 32 G — Stage-3 pre-flight gate (arm=$ARM)"
apptainer exec --bind "$SHARED" "$SIF" bash -c '
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate scenicplus
  python '"$ROOT"'/scripts/preflight.py'
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== PRE-FLIGHT PASS — now: sbatch submit_crossscore.sh ===="
else
  echo "==== PRE-FLIGHT FAILED (exit $rc) — do NOT submit the cross-scoring ===="
  echo "     Read $ROOT/out/s3_preflight_report_${ARM}.csv before deciding anything."
fi
exit "$rc"
