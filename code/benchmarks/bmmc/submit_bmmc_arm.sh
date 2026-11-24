#!/bin/bash
#SBATCH --job-name=bmmc_arm
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=long
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --propagate=NONE
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

BUNDLE="${AMLGRN_PROJECT}/benchmarks/bmmc"

[ -d "$BUNDLE" ] || { echo "MISSING BUNDLE: $BUNDLE"; exit 2; }
[ -s "$BUNDLE/submit_bmmc.sh" ] || { echo "MISSING: $BUNDLE/submit_bmmc.sh"; exit 2; }
[ "$#" -ge 1 ] || { echo "usage: sbatch submit_bmmc_arm.sh <arm> [arm ...]"; exit 2; }

echo ">> submit_bmmc_arm.sh: GRN-half resources, handing over to submit_bmmc.sh"
echo ">>   arms      : $*"
echo ">>   partition : ${SLURM_JOB_PARTITION:-unset}   cpus: ${SLURM_CPUS_PER_TASK:-unset}   mem: ${SLURM_MEM_PER_NODE:-unset}"

exec bash "$BUNDLE/submit_bmmc.sh" "$@"
