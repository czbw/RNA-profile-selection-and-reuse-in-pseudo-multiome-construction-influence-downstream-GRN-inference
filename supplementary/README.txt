AML-GRN supplementary materials

Supplementary Methods S1-S8: supplementary_methods.pdf
Supplementary Figures S1-S3: Supplementary_Figures/supplementary_figures.pdf

Supplementary Data 1: published link and metacell network membership and target agreement.
Supplementary Data 2: random-gene controls and PBMC/BMMC native-pairing benchmarks.
Supplementary Data 3: matched-reuse counts, per-stratum multiplicity checks, activity and variance.
Supplementary Data 4: matched-reuse network agreement and pAML candidate nomination.
Supplementary Data 5: adult expression scores, maturity analyses and functional annotation.
Supplementary Data 6: pediatric activity decomposition, CLP sensitivity and patient omission.

Each data folder includes an Excel workbook with selected tables and CSV files containing
the detailed results. The CSV tables retain all tested results within the stated selection,
including nonsignificant findings. Source matrices and full analysis reports are supplied
separately in source data/.

claim_index.csv links figures and manuscript claims to tables.
table_index.csv gives the source files and selection rules for every curated table.
column_sources.csv identifies the original column or calculation for each output column.
curation_manifest.json records the same table metadata and checksums.
manuscript_callout_map.csv links numbered citations to their data folders.

Oriented, centered module scores are primary. Unoriented centered scores are secondary.
The pooled cell-type display is raw, uncentered and unoriented. The Figure 5D sensitivity
removes maturity markers from regulon targets and keeps the maturity axes intact. The
earlier axis-target-exclusion analysis is separately labeled in maturity_axes_legacy.

Native-pairing networks are computational references. TF-gene correlation signs do not
establish activation or repression. Within-cell-type activity change does not imply a change
in regulatory edges. Empty cells denote missing or uncomputed values, not zero.

The tables were curated from preserved outputs without rerunning inference or statistical
tests. The code package includes scripts/curate_supplementary.py to repeat the CSV assembly.
