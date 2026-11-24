#!/usr/bin/env bash
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --job-name=pbmc10k_pairing
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=20
#SBATCH --mem=600G
#SBATCH --time=3-00:00:00
#SBATCH --output=pbmc10k_%j.out
#SBATCH --error=pbmc10k_%j.err
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
BUNDLE=$ROOT/benchmarks/pbmc
SIF=$ROOT/miniforge.sif
DS=10k

export TMPDIR=$ROOT/pbtmp
export JOBLIB_TEMP_FOLDER=$TMPDIR         # not /dev/shm: joblib memmaps there and it counts
mkdir -p "$TMPDIR" "$BUNDLE/logs"         # against the job's memory limit
[ ${#TMPDIR} -lt 100 ] || { echo "FATAL: TMPDIR too long (${#TMPDIR}b)" >&2; exit 3; }

export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"
export APPTAINERENV_AMLGRN_ROOT="$BUNDLE"
export APPTAINERENV_AMLGRN_CONFIG="config.yaml"   # common.py reads this from scripts/
export APPTAINERENV_AMLGRN_DB="$ROOT/data/reference/motifs"        # AMLGRN_ROOT is the bundle here, but the
                                                       # motif resources are shared across workflows.
export APPTAINERENV_TMPDIR="$TMPDIR"
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$TMPDIR"

echo "=================================================================="
echo " PBMC 10k pairing   job ${SLURM_JOB_ID:-interactive}   $(date)"
echo " cpus ${SLURM_CPUS_PER_TASK:-20} | mem 600G | TMPDIR $TMPDIR"
echo "=================================================================="

for p in "$SIF" "$ROOT/envs/scenicplus" "$ROOT/envs/macs3env" \
         "$BUNDLE/run_full.sh" "$BUNDLE/scripts/config.yaml" \
         "$BUNDLE/scripts/make_pairings.py" "$BUNDLE/scripts/atac_chain.py" \
         "$BUNDLE/scripts/region_sets.py" "$BUNDLE/scripts/run_grn.sh" \
         "$BUNDLE/scripts/common.py" "$BUNDLE/work/link/link_${DS}.csv"; do
  [ -e "$p" ] || { echo "FATAL: missing $p" >&2; exit 3; }
done
for f in "$BUNDLE/data/pbmc_granulocyte_sorted_${DS}_filtered_feature_bc_matrix.h5" \
         "$BUNDLE/data/pbmc_granulocyte_sorted_${DS}_atac_fragments.tsv.gz" \
         "$BUNDLE/data/pbmc_granulocyte_sorted_${DS}_atac_fragments.tsv.gz.tbi"; do
  [ -s "$f" ] || { echo "FATAL: missing/empty $f" >&2; exit 3; }
done
[ -d "$ROOT/data/reference/motifs" ] || { echo "FATAL: cisTarget feathers not found at $ROOT/data/reference/motifs" >&2; exit 3; }

echo ">> launching run_full.sh inside $SIF"
apptainer exec --bind "$ROOT" "$SIF" bash "$BUNDLE/run_full.sh" "$DS"
rc=$?
echo ">> container exited with $rc"
exit $rc
