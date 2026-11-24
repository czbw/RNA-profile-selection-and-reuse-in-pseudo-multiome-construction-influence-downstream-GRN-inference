#!/usr/bin/env bash
#SBATCH -J s3r_verify
#SBATCH --partition=bioinformatics
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
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
[ -s "$ROOT/scripts/verify.py" ] || { echo "ERROR: verify.py not uploaded."; exit 1; }

miss=0
for f in "$ROOT/out/eRegulon_AUC_Diagnosis_${ARM}_defRelapse.h5mu" \
         "$ROOT/out/eRegulon_AUC_Relapse_${ARM}_defDiagnosis.h5mu" \
         "$S2/eRegulon_AUC_Diagnosis_${ARM}.h5mu" "$S2/eRegulon_AUC_Relapse_${ARM}.h5mu" \
         "$S2/eRegulons_Diagnosis_${ARM}.tsv" "$S2/eRegulons_Relapse_${ARM}.tsv" \
         "$S2/cell_data_Diagnosis.csv" "$S2/cell_data_Relapse.csv"; do
  if [ -s "$f" ]; then printf '  OK      %-52s %s\n' "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else printf '  MISSING %s\n' "$f"; miss=1; fi
done
[ "$miss" -eq 0 ] || { echo "ERROR: cannot verify — inputs above are missing."; exit 2; }

export APPTAINERENV_CONDA_ENVS_DIRS="$SHARED/envs"
export APPTAINERENV_AMLGRN_S2="$S2"
export APPTAINERENV_AMLGRN_S3_OUT="$ROOT/out"
export APPTAINERENV_AMLGRN_S3_PAIRING="$ARM"
export APPTAINERENV_MPLCONFIGDIR="$ROOT/out/.mpl"

echo ">> bioinformatics: 4 CPUs / 64 G — Stage-3 cross-score verification (arm=$ARM)"
apptainer exec --bind "$SHARED" "$SIF" bash -c '
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate scenicplus
  python '"$ROOT"'/scripts/verify.py'
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== VERIFICATION PASSED — safe to download ===="
  ls -lh "$ROOT"/out/eRegulon_AUC_*_${ARM}_def*.h5mu 2>/dev/null
else
  echo "==== VERIFICATION FAILED (exit $rc) — do NOT download, do NOT run s3_10_prep ===="
  echo "     Read $ROOT/out/s3r_verify_${ARM}.json and the log above."
fi
exit "$rc"
