"""Create working and output directories without modifying existing files."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

bundles = [
    'paml/healthy', 'paml/malignant', 'paml/corrected', 'paml/healthy_corrected',
    'benchmarks/bmmc', 'benchmarks/pbmc', 'benchmarks/repeats',
    'scoring/published', 'scoring/corrected', 'scoring/healthy',
    'longitudinal/published', 'longitudinal/corrected',
    'longitudinal/healthy_published', 'longitudinal/healthy_corrected',
    'adult/published', 'adult/corrected', 'analysis/network_comparison',
    'preprocessing', 'modules',
]
directories = [f'{bundle}/{kind}' for bundle in bundles for kind in ('out', 'work', 'logs')]
directories += [
    'logs', 'analysis/out', 'analysis/out/gene_detection',
    'analysis/out/stratum_activity', 'analysis/out/go_stage_comparison',
    'analysis/out/marker_exclusion', 'analysis/out/patient_omission',
    'analysis/out/seven_celltypes/external',
]

for relative in directories:
    target = (ROOT / relative).resolve()
    if not target.is_relative_to(ROOT):
        raise ValueError(f'Output path is outside the repository: {relative}')
    target.mkdir(parents=True, exist_ok=True)

print(f'Prepared {len(directories)} output and working directories under {ROOT}')
