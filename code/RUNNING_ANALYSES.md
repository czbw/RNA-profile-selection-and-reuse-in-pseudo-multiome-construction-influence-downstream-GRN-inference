# Running the analyses

Run commands from the repository root unless another directory is specified. Use a writable copy, as many scripts overwrite their `out/` files. The commands retain the settings used in the study, including arguments that differ from script defaults.

## Inputs and software

To rebuild figures, run `python scripts/stage_source_data.py --source "../source data" --code-root .` from the code root. Use the supplied `figures/data/manifest.json` and `modules/out/` definitions. The manifest contains historical entries that are not used by the current figures.

To rerun the analyses, obtain the following inputs through the listed accessions and data-access procedures:

| Input | Location and required files |
|---|---|
| pAML RNA and ATAC | GSE235063 and GSE235308; metadata, indexed fragments, malignant RNA objects and published linkage tables in the paths named by the AML bundle configurations |
| PBMC multiome | Deposited 10x feature matrix and indexed fragments; `benchmarks/pbmc/work/10k/` annotated RNA, ArchR linkage/CCA files, shared ATAC chain, assignment CSVs and SCENIC+ outputs |
| BMMC multiome | GSE194122 matrix, annotations and consensus peaks; `benchmarks/bmmc/data/multiome_BMMC.h5ad`, `out/b00_cell_metadata.csv`, assignment CSVs, shared ATAC/motif outputs and per-arm networks |
| Adult expression | `adult/published/GSE327386_raw_count_matrix.txt.gz`, original gene identifiers, and `data/reference/lambo_cell_markers.xlsx` for the maturity-axis definitions |
| Motif and genome resources | GRCh38 cisTarget ranking databases, motif-to-TF annotations, genome annotation, chromosome-size files and TF lists specified in each bundle configuration |
| Functional annotation | Files under `analysis/ref/`, including the Lambert TF table and the GO/Reactome/Hallmark GMT libraries; use the versions specified for the study |
| Cross-stage activity | Malignant and healthy AUCell matrices for both stage vocabularies, with patient/stage/cell-type metadata and the same composite cell identifiers as the original analysis |

The source-data folder does not include all raw data, motif databases or intermediate matrices required for network inference.

Python analyses use NumPy, pandas, SciPy, statsmodels, scikit-learn, anndata, h5py, PyYAML, ctxcore, networkx and openpyxl. Chromatin processing and network inference require the specified SCENIC+, pycisTopic, pycisTarget and MACS3 environments. R figures use ggplot2, patchwork, jsonlite, scales, grid and ragg. The graphical abstract requires ReportLab, Pillow, Arial fonts and Poppler. The supplied YAML files specify dependencies but do not fully record the software installed for each historical run.

Before submitting cluster jobs, set the allocation, partition, email, container and resource paths. From the repository root in bash, run `export AMLGRN_PROJECT="$PWD"` before invoking `sbatch`; launchers use this location for shared resources. Motif resources are stored in `data/reference/motifs/`, and workflow-specific configurations may require links to that directory. Run shell scripts with bash and preserve LF line endings. Bash line continuations cannot be pasted directly into PowerShell.

## Network and activity generation

1. Stage pAML metadata and fragments using `preprocessing/scripts/extract_atac.sh`, `audit_atac.py` and `prepare_fragments.sh`. Shared resources are downloaded by `download_databases.sh`. The full healthy and malignant bundles consume these staged resources; paths must match their configurations.
2. Use `paml/healthy/` for the shared healthy construction and `paml/malignant/` for malignant link/metacell networks. Prepare malignant RNA with `paml/malignant/scripts/build_malignant_h5ad.py`. Run the bundle driver with a specific stage, for example `bash run_all.sh Diagnosis` from the relevant bundle on the configured cluster. Use a separate invocation for `Relapse`. The two stages can also be run together by omitting the stage argument.
3. For PBMC, run preparation/annotation, ArchR linking and pairing construction through `prepare.py`, `annotate.py`, `archr_link.R` and `make_pairings.py`. The configured `run_full.sh` builds one shared chromatin chain and runs selected arms. `build_arm_rna.py` stages alternative assignments; `run_grn.sh` handles each selected arm against the same chain.
4. For BMMC, `census.py`, `export.py`, `export_gene_activity.py` and the `infer_pairings.py`, `infer_stratified_pairings.py` and `infer_capped_pairings.py` prepare inputs and assignment tables. `arms.py` makes random and reuse-matched comparators. `run_bmmc.sh` and `run_grn_bmmc.sh` share a single chromatin chain across arms. Do not independently rebuild each arm's chromatin input.
5. Build corrected malignant and healthy assignments with `paml/corrected/` and `paml/healthy_corrected/`. The corrected malignant census requires the `--malignant-rna-only` option below. Corrected malignant network inference uses `run_grn.sh`; healthy corrected profiles are scored against the malignant target vocabularies.
6. Cross-stage rescoring uses `scoring/published/`, `scoring/corrected/`, `scoring/healthy/` and `paml/healthy_corrected/scripts/crossscore.sh`. Their checks verify cell identities and target-vocabulary correspondence.

## Benchmark commands

Prepare the input tables and matrices before running these commands. Replace `<project-root>` and output-directory placeholders with local paths.

| Analysis | Command and working directory |
|---|---|
| Eligible malignant RNA census | `python paml/corrected/scripts/feasibility_census.py --root <project-root> --malignant-rna-only` |
| Corrected malignant assignments | `python paml/corrected/scripts/capacity_assign.py --root <project-root> --seed 555` |
| BMMC unrestricted inference | `python benchmarks/bmmc/scripts/infer_pairings.py --knn 100` |
| BMMC matched-reuse checks (1.406) | `python benchmarks/bmmc/scripts/verify_matched_reuse.py` |
| PBMC unrestricted one-to-one inference | `python benchmarks/pbmc/scripts/link_1to1_assign.py --dataset 10k --knn 240` |
| PBMC within-subcluster one-to-one inference | `python benchmarks/pbmc/scripts/link_subcluster_1to1_assign.py --dataset 10k --dense-check` |
| PBMC full random assignment report | `python benchmarks/pbmc/scripts/make_random_arms.py --dataset 10k --seed 555 --min-stratum 0` |
| PBMC restricted random assignment report | `python benchmarks/pbmc/scripts/make_random_arms.py --dataset 10k --seed 555` |
| PBMC activity agreement | `python benchmarks/pbmc/scripts/activity_agreement.py --dataset 10k --seeds 555` |
| PBMC group-mean profiles | `python benchmarks/pbmc/scripts/group_mean_activity.py --dataset 10k --seeds 555` |
| BMMC activity agreement | From `benchmarks/bmmc/scripts/`: `python activity_agreement.py --out-dir ../out --metadata ../out/b00_cell_metadata.csv` |
| BMMC group-mean profiles | From `benchmarks/bmmc/scripts/`: `python group_mean_activity.py` |
| Matched-reuse donor-stratum mean activity | `python analysis/activity/stratum_mean_activity.py` |
| RNA nearest-neighbor comparator | `python analysis/activity/nearest_neighbor_comparator.py` |
| Identical-input replicate comparison | `python benchmarks/repeats/scripts/compare_replicates.py --src benchmarks/pbmc/out --rep benchmarks/repeats/out` |

The same script summarizes PBMC and BMMC network agreement:

```text
python benchmarks/bmmc/scripts/network_agreement.py --out-dir benchmarks/pbmc/out --stage pbmc10k --arms link link_1to1 link_subcluster_1to1 pseudomultiome pseudomultiome_k10 rand_subcluster_matched rand_subcluster_1to1 rand_global --unpaired pseudomultiome pseudomultiome_k10 --report-name trackb_rungs_pbmc10k_all.json
python benchmarks/bmmc/scripts/network_agreement.py --out-dir benchmarks/bmmc/out --stage bmmc --arms rand_celltype_1to1 rand_celltype_x_donor_1to1 rand_global link link_1to1 link_x_celltype_donor_1to1 link_x_celltype_donor rand_celltype_donor_matched link_x_celltype_donor_cap14 rand_celltype_donor_matched_cap14 pseudomultiome pseudomultiome_k10 --unpaired pseudomultiome pseudomultiome_k10 --report-name trackb_rungs_bmmc_all.json
```

The capacity sweep was produced by two runs:

```text
python benchmarks/bmmc/scripts/infer_capped_pairings.py --sweep 1.0 --outdir <dirA> --space-cache <cache>
python benchmarks/bmmc/scripts/infer_capped_pairings.py --sweep 0.2,0.3,0.4,0.5,0.6,0.7,0.8 --outdir <dirB> --space-cache <cache>
```

The final summary was assembled manually by appending the first run's result to the second run's entries and retaining the second run's other top-level fields. The `sweep_note` field was added by hand. No automated merge script is available.

## pAML summaries and module definition

Generate the Figure 1 summaries with `analysis/pairing/pairing_target_sets.py`, `celltype_profiles.py`, `celltype_labels.py`, `shared_regulon_activity.py`, `cell_count_sensitivity.py` and `activity_target_agreement.py`. The random-gene-set control is `random_gene_sets.py`. `tf_families.py` and `pairing_enrichment.py` provide TF-family assignments and expression-based gene backgrounds for subsequent analyses.

For each pAML construction, `prep.py`/`prep.py` produce patient-cell-type pseudobulks and paired-unit tables. The model, offset, recalibration and candidate steps are `model.py`, `offset.py`, `recalibrate.py` and `candidates.py`; the controls are in `controls.py`. Complete healthy-cell preparation, modeling and selectivity analyses before nomination, and the within-cell-type decomposition step `decomposition.py` before `nomination.py`. Corrected malignant analyses use the scripts under `longitudinal/corrected/scripts/`. Healthy corrected selectivity is `longitudinal/healthy_corrected/scripts/selectivity.py`.

`analysis/network_comparison/scripts/compare_networks.py` generates the regulon catalogue required for nomination. `modules/scripts/module_definitions.py --mode build` writes `corrected_module.json` and related definitions. Use the 55-member T3_repair module for the reported analyses. The older 34-member `adult/corrected/out/module_definition.json` supplies background targets and orientations to earlier processing steps; it must not replace the 55-member definition.

Some scripts require inputs for additional calculations that are not reported in the manuscript. The published link nomination script reads existing K562 tables for historical output fields, although these fields do not determine the reported membership. Supply those tables when running the script. The ribosomal-control script reads malignant gene-detection tables from `analysis/out/gene_detection/`, produced by `analysis/activity/gene_detection.py`. This calculation uses RNA expression and does not perform ENCODE validation. The enrichment background is defined separately by `pairing_enrichment.py`.

## Adult scores and sensitivity analyses

`adult/published/scripts/feasibility.py` checks the adult cohort. The target and background definitions require `adult/published/scripts/build_module.py`, `adult/corrected/scripts/build_module.py` and `adult/corrected/scripts/prepare_candidate_table.py`. `adult/corrected/scripts/score.py` calculates AUCell, singscore-based and mean-standardized-expression scores over the full corrected-arm background. The following tests use the 55-member module:

```text
python adult/shared/coverage.py --module modules/out/corrected_module.json --arm repair --tag T3_repair
python adult/shared/nullB_pool.py --module modules/out/corrected_module.json --arm repair --tag T3_repair
python adult/shared/test.py --module modules/out/corrected_module.json --arm repair --tag T3_repair
python adult/shared/robustness.py --module modules/out/corrected_module.json --arm repair --tag T3_repair
python adult/shared/axis_lambo.py --module modules/out/corrected_module.json --arm repair --tag T3_repair
python analysis/sensitivity/exclude_maturity_markers.py
python analysis/sensitivity/verify_marker_exclusion.py
```

The first five commands specify the module, construction and output tag. Marker exclusion removes shared maturity-marker genes from the scored regulon targets while retaining the intact maturity axes.

Module characterization uses `regulon_characterization.py`, `module_characterization.py`, `module_network.py`, `ribosomal_controls.py`, `patient_unit_adjustment.py`, `module_decomposition.py` and `module_figure_data.py`. The stagewise annotation comparison is `go_stage_comparison.py`, followed by `summarize_go_stage_comparison.py`. These scripts support Figure 5 and Supplementary Data 5-6.

Seven-cell-type sensitivity uses `seven_celltype_nomination.py --repo <project-root> --out <output-directory>`, followed by `seven_celltype_module.py` with the same arguments. Score the resulting module with the adult-test scripts, specifying its output path and the sensitivity tag. Patient omission uses `patient_omission.py`; its `report_patient_omission.py` generates the S1/S2 source plots from `analysis/out/patient_omission/`.

## Reproducibility limits

- The original script that produced `trackA_L6_10k_subcluster_arm.json` was not saved. The result file is supplied; the later reconstructed script is excluded.
- The matched-random report at reuse 1.406 (`cap14`) was produced by `arms.py` and subsequently renamed. The original rename command was not recorded.
- The BMMC group-mean script now scores eleven constructions, whereas the supplied result contains six. Earlier checks found identical numerical values for the six shared constructions; the complete files therefore differ.
- Exact script versions and installed software environments are not fully recorded for each historical cluster job.

## Output folders and source logs

Run `python scripts/prepare_outputs.py` before the analyses. It creates output and working directories; raw data must be supplied separately. The code archive omits empty directories.

`analysis/pairing/aml_metacell_counts.py` reads the original SCENIC+ job logs from `data/logs/diagnosis_metacells.log` and `data/logs/relapse_metacells.log`. The supplied figure inputs already contain the extracted counts.

Internal result identifiers such as `link`, `repair`, `T3_repair` and `cap14` are retained in arguments and data fields to match the recorded results. They denote the published link construction, corrected construction, 55-regulon corrected module and matched-reuse level near 1.406, respectively. They do not identify separate code versions.
