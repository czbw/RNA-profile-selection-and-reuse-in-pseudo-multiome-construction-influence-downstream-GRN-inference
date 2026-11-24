#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v conda >/dev/null || { echo "ERROR: conda not found. module load / source your conda first."; exit 1; }
# mamba is much faster if available; conda works too.
SOLVER=conda; command -v mamba >/dev/null && SOLVER=mamba

echo ">> creating env 'scenicplus' (scenicplus 1.0a2 / pycisTopic 2.0a / pyscenic 0.12.1 + pins)"
$SOLVER env create -f "$HERE/environment.yml"

echo ">> creating env 'macs3env' (separate environment)"
$SOLVER env create -f "$HERE/environment.macs3.yml"

echo ">> verifying scenicplus imports"
conda run -n scenicplus python -c "import scenicplus, pycisTopic, pyscenic, mudata, polars; print('scenicplus', scenicplus.__version__)"
echo ">> macs3 at: $(conda run -n macs3env which macs3)"
echo ">> DONE. Both envs built. Next: edit run_all.sh #SBATCH headers, then sbatch run_all.sh"
