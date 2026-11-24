Supplementary Data 2. Activity controls and native-pairing benchmarks

CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.
The table index identifies source files and selection rules. Source files are preserved under source data/.
source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.
No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.

Tables:
random_gene_controls.csv - Expression- and size-matched random gene controls for Figure 1.
random_control_summary.csv - Summary of the random gene controls cited in Results.
pbmc_network_agreement.csv - PBMC native-reference network comparisons (Figure 2A,B).
pbmc_network_details.csv - Counts and additional network metrics for PBMC.
bmmc_network_agreement.csv - BMMC native-reference network comparisons (Figure 2A,B).
bmmc_network_details.csv - Counts and additional network metrics for BMMC.
pbmc_cell_activity.csv - Per-cell activity agreement (Figure 2C).
pbmc_stratified_cell_activity.csv - Per-cell activity agreement (Figure 2C).
bmmc_cell_activity.csv - Per-cell activity agreement (Figure 2C).
bmmc_assignment_accuracy.csv - Assignment accuracy and label agreement in BMMC.
rna_neighbor_comparator.csv - RNA nearest-neighbor comparator (Figure 2C).
bmmc_cell_type_activity.csv - Underlying cell-type or regulon activity comparisons.
bmmc_regulon_activity.csv - Underlying cell-type or regulon activity comparisons.
bmmc_group_mean_activity.csv - Underlying cell-type or regulon activity comparisons.

Definitions

The random-gene control uses N1 from the preserved report: gene sets matched in size
and mean-expression-bin counts, excluding the regulon's own targets. N2 and N3 are
additional controls retained in source data, not the control quoted in the Results.
L1 compares top region-gene pairs and their correlation coefficients. L2 compares top
TF-gene edges. L3 compares shared-regulon targets; L4 compares signed eRegulon identities;
L5 compares TF identities. Native-pairing networks are inferred computational references.
Arm identifiers match the source reports: link denotes inferred pairing in these benchmarks,
1to1 denotes one-to-one assignment, and subcluster or celltype_donor denotes restriction
to the stated strata. pseudomultiome denotes metacell aggregation. These benchmark link
identifiers do not refer to the pAML published link construction.
