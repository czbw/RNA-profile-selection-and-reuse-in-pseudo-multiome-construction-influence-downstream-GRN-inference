# Code organization

Scripts are grouped by dataset and analysis. Paths, imports, configuration defaults and run instructions use the same layout. Empty output directories are omitted; `scripts/prepare_outputs.py` creates the required directories.

The gene-detection script is under `analysis/activity/`. Adult scoring utilities shared by both constructions are under `adult/shared/`. Preprocessing inputs and derived files are stored outside `preprocessing/scripts/`.

The pAML drivers no longer call the absent `compare.py` development script. The PBMC preflight and staging scripts refer to the supplied runners. Local defaults are derived from the repository location; cluster jobs still require site-specific resource and environment settings.

Metacell counts are extracted from two unchanged SCENIC+ job logs supplied in the figure-input archive. Their per-cell-type counts match the recorded results. Numerical analysis settings, result identifiers and supplied result tables are unchanged.

Figures 1-5, Supplementary Figures S1-S3 and the graphical abstract were rebuilt and matched the current submission visuals pixel for pixel. Marker-exclusion checks, Python/R/Bash syntax checks and workflow-root checks passed. Full network inference and cluster execution were not repeated.

## Source outputs and supplementary tables

Preserved outputs are supplied in the separate `source data` folder. Its manifest replaces
the figure-input manifest. `scripts/stage_source_data.py` copies these outputs into their
expected working paths. `scripts/curate_supplementary.py` assembles the numbered supplementary
CSV tables from the same files. These are packaging utilities, not missing analysis producers.
