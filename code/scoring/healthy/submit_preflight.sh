#!/usr/bin/env bash
#SBATCH -J s3b_preflight
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

ROOT="${AMLGRN_PROJECT}/scoring/healthy"
SHARED="${AMLGRN_PROJECT}"
HEALTHY="${AMLGRN_HEALTHY:-$SHARED/out}"
MALIG="$SHARED/paml/malignant/out"
SIF="$SHARED/miniforge.sif"
mkdir -p "$ROOT/logs" "$ROOT/out"

[ -f "$SIF" ]                    || { echo "ERROR: $SIF not found."; exit 1; }
[ -d "$SHARED/envs/scenicplus" ] || { echo "ERROR: $SHARED/envs/scenicplus missing."; exit 1; }
[ -s "$ROOT/scripts/preflight.py" ] || { echo "ERROR: preflight.py not uploaded."; exit 1; }

# The healthy mudata are Stage-1 outputs; if run_grn.sh wrote them somewhere other than out/,
# report that explicitly rather than failing later with a less clear message.
miss=0
for f in "$HEALTHY/mudata_Diagnosis_link.h5mu" "$HEALTHY/mudata_Relapse_link.h5mu" \
         "$MALIG/mudata_Diagnosis_link.h5mu"   "$MALIG/mudata_Relapse_link.h5mu" \
         "$MALIG/eRegulons_Diagnosis_link.tsv" "$MALIG/eRegulons_Relapse_link.tsv"; do
  if [ -s "$f" ]; then printf '  OK      %-52s %s\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else printf '  MISSING %s\n' "$f"; miss=1; fi
done
if [ "$miss" -ne 0 ]; then
  echo
  echo "ERROR: inputs missing. If the HEALTHY mudata are not in $HEALTHY, find them with:"
  echo "         find $SHARED -maxdepth 3 -name 'mudata_*_link.h5mu' -ls"
  echo "       then re-run with:  AMLGRN_HEALTHY=<dir> sbatch submit_preflight.sh"
  exit 2
fi

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
export APPTAINERENV_AMLGRN_HEALTHY="${AMLGRN_HEALTHY:-$HEALTHY}"
export APPTAINERENV_AMLGRN_MALIG="$MALIG"
export APPTAINERENV_AMLGRN_S3B_OUT="$ROOT/out"

echo ">> bioinformatics: 2 CPUs / 32 G — Stage-3B pre-flight gate"
apptainer exec --bind "$SHARED" "$SIF" bash -c '
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate scenicplus
  cd '"$ROOT"'/scripts
  python preflight.py'
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== PRE-FLIGHT PASS ===="
  echo "     Read the MEMORY GUIDANCE lines above, set --mem in submit_crossscore.sh to match,"
  echo "     then: sbatch submit_crossscore.sh"
else
  echo "==== PRE-FLIGHT FAILED (exit $rc) — do NOT submit the cross-scoring ===="
  echo "     Read $ROOT/out/s3b_preflight_report.csv before deciding anything."
fi
exit "$rc"
