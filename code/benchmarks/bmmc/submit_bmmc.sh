#!/bin/bash
#SBATCH --job-name=bmmc_ladder
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=40
#SBATCH --mem=2400G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --propagate=NONE
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_PROJECT}"
BUNDLE=$ROOT/benchmarks/bmmc
SIF=$ROOT/miniforge.sif

export TMPDIR=$ROOT/bmtmp
mkdir -p "$TMPDIR" "$BUNDLE/logs" "$BUNDLE/out" "$BUNDLE/work"

export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"          # else EnvironmentNameNotFound
export APPTAINERENV_AMLGRN_ROOT="$BUNDLE"
export APPTAINERENV_AMLGRN_STAGE="bmmc"
export APPTAINERENV_AMLGRN_CONFIG="config.yaml"
export APPTAINERENV_AMLGRN_DB="$ROOT/data/reference/motifs"           # GRCh38 feathers, reused from the AML run
export APPTAINERENV_AMLGRN_NCPU="${SLURM_CPUS_PER_TASK:-20}"
export APPTAINERENV_AMLGRN_SEED="${AMLGRN_SEED:-555}"
export APPTAINERENV_AMLGRN_ARMS="${*:-truth rand_celltype_1to1 rand_celltype_x_donor_1to1 rand_global}"
export APPTAINERENV_AMLGRN_STOP_AFTER="${AMLGRN_STOP_AFTER:-}"
export APPTAINERENV_TMPDIR="$TMPDIR"
export APPTAINERENV_JOBLIB_TEMP_FOLDER="$TMPDIR"

# PREFLIGHT. Checked here, on the host, before the container is entered, so a missing file fails
# instantly with its name instead of as a ModuleNotFoundError after the scheduler has granted the
# allocation.
for p in "$SIF" "$ROOT/envs/scenicplus" \
         "$BUNDLE/run_bmmc.sh" \
         "$BUNDLE/scripts/run_grn_bmmc.sh" \
         "$BUNDLE/scripts/chain_guard.sh" \
         "$BUNDLE/scripts/atac_chain.py" \
         "$BUNDLE/scripts/region_sets.py" \
         "$BUNDLE/scripts/build_arm_rna.py" \
         "$BUNDLE/scripts/config.yaml"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# the source matrix and the barcode-remap table
for p in "$BUNDLE/data/multiome_BMMC.h5ad" "$BUNDLE/uploads/b00_cell_metadata.csv"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# the GRCh38 cisTarget databases — reused, never staged for this bundle
for p in "$ROOT/data/reference/motifs/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather" \
         "$ROOT/data/reference/motifs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
         "$ROOT/data/reference/motifs/hg38_genome_annotation.tsv" \
         "$ROOT/data/reference/motifs/hg38_chromsizes_3col.tsv"; do
  [ -e "$p" ] || { echo "MISSING: $p"; exit 2; }
done

# every requested arm that needs an uploaded assignment table must have one. `truth` is barcode
# identity and the metacell arms have no per-cell assignment, so both are exempt by design.
SEED_USED="${AMLGRN_SEED:-555}"
for ARM in ${*:-truth rand_celltype_1to1 rand_celltype_x_donor_1to1 rand_global}; do
  case "$ARM" in
    truth|pseudomultiome|pseudomultiome_k*) continue ;;
  esac
  p="$BUNDLE/uploads/arm_${ARM}_bmmc_seed${SEED_USED}_full.csv"
  [ -e "$p" ] || { echo "MISSING assignment table: $p"; exit 2; }
done

# A CRLF in any of these scripts makes bash die on `set -o pipefail` with an opaque error.
# Cheap to check, so check it rather than trust it.
for p in "$BUNDLE/run_bmmc.sh" "$BUNDLE/submit_bmmc.sh" "$BUNDLE/scripts/run_grn_bmmc.sh" \
         "$BUNDLE/scripts/chain_guard.sh"; do
  if LC_ALL=C grep -q $'\r' "$p" 2>/dev/null; then
    echo "CRLF DETECTED in $p -- run: sed -i 's/\r$//' $p"; exit 2
  fi
done

echo ">> preflight OK: scripts, databases, source matrix and assignment tables all present"
echo ">> entering $SIF once"
echo ">>   arms : ${*:-truth rand_celltype_1to1 rand_celltype_x_donor_1to1 rand_global}"
echo ">>   seed : $SEED_USED   stop_after: '${AMLGRN_STOP_AFTER:-}'"
apptainer exec --bind "$ROOT" "$SIF" bash "$BUNDLE/run_bmmc.sh"
