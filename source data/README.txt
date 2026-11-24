AML-GRN source data

This folder preserves analysis outputs and the inputs required by the figure builders.
File contents have not been edited. manifest.csv records their paths and SHA-256 checksums.
manifest.json also records original locations and their former supplementary-data locations.
The numbered Supplementary Data packages contain selected, organized tables drawn from these files.

The source files include full result tables, score matrices, random-draw arrays and job logs.
They do not include all raw sequencing data, motif resources or intermediate matrices needed
to rerun network inference. Dataset accessions and analysis inputs are documented in the code.

Some source reports contain earlier analysis labels or superseded interpretations. Use the
manuscript and curated supplementary tables to identify the reported analyses. In particular,
the oriented centered module score is primary, and Figure 5D retains intact maturity axes
while removing overlapping markers from regulon targets. The capacity-sweep report was
already a manual merge of two runs; it is identified as such in the manifest.

To stage these files for the figure scripts, run from the code folder:
python scripts/stage_source_data.py --source "../source data" --code-root .

The helper verifies checksums and copies each file to its recorded repository path. It refuses
to replace a different existing output. Run analyses in a working copy so that this source
archive remains unchanged. No code is executed from this folder.
