# AML-GRN

Code for evaluating RNA-ATAC pairing constructions, comparing inferred gene regulatory networks, and analysing diagnosis-to-relapse activity changes in AML.

## Start here

To reproduce the figures, place the supplied `source data` folder beside this code folder. Then run:

```sh
python scripts/stage_source_data.py --source "../source data" --code-root .
python scripts/prepare_outputs.py
Rscript figures/main/make_fig01_construction.R
Rscript figures/main/make_fig02_ladder.R
Rscript figures/main/make_fig03_concentration.R
Rscript figures/main/make_fig04_consequence.R
Rscript figures/main/make_fig05_module.R
Rscript figures/supplementary/make_figS_reuse_native_agreement.R
python analysis/sensitivity/report_patient_omission.py
python figures/graphical_abstract/build_graphical_abstract.py
```

The patient-omission reporter produces Supplementary Figures S1 and S2. The other builders write their figures beside their scripts. See [figure instructions](figures/README.md) for software requirements.

For the analysis sequence, required inputs and study-specific command arguments, see [RUNNING_ANALYSES.md](RUNNING_ANALYSES.md). Network inference requires the raw data and configured SCENIC+ environments; the supplied outputs alone are insufficient.

## Repository layout

| Folder | Purpose |
| --- | --- |
| `preprocessing/` | pAML fragment preparation and motif-resource downloads |
| `paml/` | Healthy and malignant pAML processing and corrected pairing construction |
| `benchmarks/` | PBMC and BMMC pairing benchmarks and repeated-inference checks |
| `scoring/` | Cross-stage AUCell scoring with diagnosis and relapse vocabularies |
| `longitudinal/` | Patient-cell-type summaries, statistical comparisons and candidate nomination |
| `adult/` | Adult AML expression scoring, null tests and maturity-axis analyses |
| `modules/` | Module definitions and the script that constructs them |
| `analysis/` | Pairing summaries, activity comparisons, functional analyses and sensitivities |
| `figures/` | Main figures, supplementary figures and graphical abstract |
| `data/` | Input-location instructions; raw data are supplied separately |
| `scripts/` | Repository setup utilities |
| `provenance/` | Source-output links to their producing scripts |

Within `longitudinal/`, `published` and `corrected` refer to the two malignant constructions; `healthy_published` and `healthy_corrected` contain their healthy-cell comparisons. In `adult/`, utilities under `shared/` accept either construction, while `published/` and `corrected/` prepare construction-specific inputs.

Empty `out/`, `work/` and `logs/` directories are omitted from the code archive. `prepare_outputs.py` creates them when needed. The populated `modules/out/` directory contains supplied module definitions and should be retained.

Recorded data identifiers, including `repair` and `T3_repair`, remain unchanged inside result tables and command arguments. Their meaning is documented in the analysis instructions. Source locations and checksums are listed in `release_manifest.json`; source-data checksums are listed in `../source data/manifest.json`.

## Reproducibility limits

The supplied results support figure reconstruction. Upstream inference also requires raw matrices, fragments, motif resources and cluster configuration. One original PBMC activity-report producer was not saved, and the BMMC capacity-sweep summary was assembled manually from two runs. The analysis instructions describe these gaps and the available checks. No missing producer or result has been reconstructed as though it were an original study output.

## Supplementary tables

`source data/` contains preserved outputs. The separate `supplementary/` folder contains
the six curated data packages, Supplementary Methods and Supplementary Figures.
`supplementary/table_index.csv` records the source files and selection rules for each table.

To regenerate the curated CSV tables in a new folder:

```sh
python scripts/curate_supplementary.py --source "../source data" --output ../supplementary_rebuilt
```

This utility uses the Python standard library. It selects recorded results and recomputes
only descriptive counts and ratios; it does not rerun statistical tests or network inference.
The Excel workbooks are reading copies of selected CSV tables. Figures continue to read
the preserved analysis outputs, so editing a supplementary table does not alter a figure.
