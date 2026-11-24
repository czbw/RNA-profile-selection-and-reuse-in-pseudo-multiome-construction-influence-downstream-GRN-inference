#!/bin/bash
#SBATCH --job-name=aml_healthy_repair
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=10
#SBATCH --mem=600G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
# The healthy-cell run sits at the bundle root on the cluster: its chain, outputs and
# staged RNA are at $ROOT/{out,work,data/rna}. All other bundles are nested subdirectories.
ARMA=$ROOT                                # owns the HEALTHY ATAC chain + healthy out/ + data/rna
STAGE2=$ROOT/paml/malignant       # owns the malignant repair vocabularies
BUNDLE=$ROOT/paml/healthy_corrected     # owns these scripts + the repaired healthy cell_data
SIF=$ROOT/miniforge.sif

export TMPDIR=$ROOT/hrtmp
mkdir -p "$TMPDIR" "$BUNDLE/logs" "$BUNDLE/out"

STAGES="${*:-Diagnosis Relapse}"

export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"          # else EnvironmentNameNotFound
export APPTAINERENV_AMLGRN_HEALTHY_OUT="$ARMA/out"
export APPTAINERENV_AMLGRN_HEALTHY_WORK="$ARMA/work"
export APPTAINERENV_AMLGRN_HEALTHY_RNA="$ARMA/data/rna"
export APPTAINERENV_AMLGRN_MALIG="$STAGE2/out"
export APPTAINERENV_AMLGRN_HR_BUNDLE="$BUNDLE"
export APPTAINERENV_AMLGRN_HR_OUT="$BUNDLE/out"
export APPTAINERENV_AMLGRN_STAGES="$STAGES"
export APPTAINERENV_AMLGRN_DEFS="Diagnosis Relapse"
export APPTAINERENV_AMLGRN_CELL_ID_SEP="___"
export APPTAINERENV_AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-10}"
export APPTAINERENV_TMPDIR="$TMPDIR"
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$TMPDIR"
export APPTAINERENV_MPLCONFIGDIR="$BUNDLE/out/.mpl"

# Preflight: check all inputs outside the container so a missing file fails fast with a clear message.
for p in "$SIF" "$ROOT/envs/scenicplus" \
         "$BUNDLE/run_healthy.sh" \
         "$BUNDLE/scripts/crossscore.sh" \
         "$BUNDLE/scripts/build_paired_rna.py"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# the healthy RNA the builder resolves donors against: all three files, even for a two-stage run
for p in "$ARMA/data/rna/healthy_diagnosis.h5ad" \
         "$ARMA/data/rna/healthy_remission.h5ad" \
         "$ARMA/data/rna/healthy_relapse.h5ad"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# the malignant repair vocabularies this arm is scored under
for D in Diagnosis Relapse; do
  [ -e "$STAGE2/out/eRegulons_${D}_repair.tsv" ] \
    || { echo "MISSING: $STAGE2/out/eRegulons_${D}_repair.tsv"; exit 2; }
done

# the healthy ATAC chain, and the uploaded assignment
CHAIN_MISSING=0
for STAGE in $STAGES; do
  [ -e "$ARMA/out/cistopic_obj_${STAGE}.pkl" ] \
    || { echo "MISSING CHAIN: $ARMA/out/cistopic_obj_${STAGE}.pkl"; CHAIN_MISSING=1; }
  [ -e "$BUNDLE/uploads/cell_data_${STAGE}.healthy_repaired.csv" ] \
    || { echo "MISSING: $BUNDLE/uploads/cell_data_${STAGE}.healthy_repaired.csv"; exit 2; }
done
if [ "$CHAIN_MISSING" -ne 0 ]; then
  echo
  echo "The healthy cisTopic object(s) above are absent. That decides the size of this job and it"
  echo "is not something to work around here:"
  echo "  * present -> this job is hours (GEX build + prepare_GEX_ACC + 4 AUCell runs)"
  echo "  * absent  -> the healthy ATAC chain must be rebuilt first by the Arm-A run at \$ROOT"
  echo "               (atac_chain.py: pseudobulk -> macs3 -> consensus -> cisTopic -> topics),"
  echo "               which is a multi-day job on its own. Do NOT bolt it onto this script."
  echo "Check what actually survived with:"
  echo "  ls -lh $ARMA/out/cistopic_obj_*.pkl $ARMA/out/menr_cistarget_*.hdf5; ls -d $ARMA/out/region_sets_*"
  exit 3
fi
echo ">> preflight OK: scripts, healthy RNA, healthy ATAC chain, malignant repair defs, assignment tables"

echo ">> entering $SIF once; healthy stages: $STAGES"
apptainer exec --bind "$ROOT" "$SIF" bash "$BUNDLE/run_healthy.sh"
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==== HEALTHY RE-PAIRING COMPLETE ===="
  ls -lh "$BUNDLE"/out/eRegulon_AUC_H*_repair_defM*.h5mu 2>/dev/null
else
  echo "==== HEALTHY RE-PAIRING FAILED (exit $rc) ===="
fi
exit "$rc"
