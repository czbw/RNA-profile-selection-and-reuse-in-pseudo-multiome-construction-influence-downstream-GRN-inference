# Analysis inputs

Run `python scripts/stage_source_data.py --source "../source data" --code-root .` from the code root to place the supplied outputs in their expected locations. These files are sufficient for the figure builders, but not for repeating upstream network inference.

| Location | Input |
| --- | --- |
| `paml/GSE235063_RAW/` | pAML RNA source files |
| `paml/atac_meta/`, `paml/fragments_raw/`, `paml/fragments/` | pAML metadata and fragment files prepared by `preprocessing/scripts/` |
| `pbmc/` | PBMC multiome feature matrices and fragments |
| `rna_panel/` | `panel_genes.txt` and `panel_meta.csv` for gene-detection diagnostics |
| `reference/lambo_cell_markers.xlsx` | Lambo et al. Supplementary Table S3 used for maturity markers |
| `reference/motifs/` | cisTarget databases and related resources |
| `logs/diagnosis_metacells.log`, `logs/relapse_metacells.log` | Original SCENIC+ metacell-count logs supplied in the source-data folder |
| `reference_support/` | Existing K562 tables read by ancillary fields in the nomination scripts; these fields do not select the reported 55-regulon module |

The pAML ATAC source archive, `GSE235308_RAW.tar`, is expected at the repository root. Adult counts are expected at `adult/published/GSE327386_raw_count_matrix.txt.gz`. BMMC input objects and workflow-specific matrices belong in the locations specified in `RUNNING_ANALYSES.md` and the workflow configuration files.

Raw data, large databases and most intermediate matrices are not included in the code package. Obtain the source data through the stated accessions and access procedures. Keep the study's gene identifiers and annotation versions when reproducing its results.
