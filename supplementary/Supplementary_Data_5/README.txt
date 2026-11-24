Supplementary Data 5. Adult expression scores and functional annotation

CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.
The table index identifies source files and selection rules. Source files are preserved under source data/.
source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.
No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.

Tables:
adult_module_scores.csv - Oriented, background-centered adult module changes (Figure 5B), including alternative scoring methods.
adult_nulls.csv - Null distributions for the primary adult module score (Figure 5B).
adult_patient_changes.csv - Adult expression-score changes, sample identifiers and target coverage (Figure 5A,B).
adult_regulons_diagnosis.csv - Adult expression-score changes, sample identifiers and target coverage (Figure 5A,B).
adult_regulons_relapse.csv - Adult expression-score changes, sample identifiers and target coverage (Figure 5A,B).
adult_target_coverage.csv - Adult expression-score changes, sample identifiers and target coverage (Figure 5A,B).
adult_samples.csv - Adult expression-score changes, sample identifiers and target coverage (Figure 5A,B).
module_membership.csv - Module membership and pediatric orientations within the 201-regulon background.
module_targets.csv - Target definitions for the module and its centering background.
maturity_marker_exclusion\summary.csv - Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).
maturity_marker_exclusion\target_retention.csv - Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).
maturity_marker_exclusion\axis_genes.csv - Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).
maturity_marker_exclusion\per_patient_changes.csv - Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).
maturity_marker_exclusion\per_regulon_target_audit.csv - Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).
maturity_marker_exclusion\module_change_comparison.csv - Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).
maturity_axis_changes.csv - Diagnosis-to-relapse change in the intact axes. Larger scores indicate primitive marker expression.
maturity_axes_legacy\axis_target_exclusion.csv - Earlier sensitivity excluding module targets from the axes, retained for Supplementary Methods S7 only.
maturity_axis_random_controls.csv - Size- and detection-matched random-axis controls (Supplementary Methods S7).
pooled_cell_type_activity.csv - Secondary pooled raw, uncentered and unoriented activity and six-cell-type shares (Figure 5C).
signed_target_enrichment.csv - All tested terms for pooled signed target sets, before and after ribosomal depletion (Figure 5E).
tf_signed_target_counts.csv - Target counts after separating each TF by its TF-gene correlation sign (Figure 5F).
tf_signed_enrichment.csv - All tested terms for signed TF target sets (Figure 5F).
ribosomal_genes.csv - Ribosomal protein genes excluded in the reported depletion sensitivity.
stagewise_function\enrichment_counts.csv - Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).
stagewise_function\gene_membership.csv - Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).
stagewise_function\go_enrichment_all_terms.csv - Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).
stagewise_function\go_exclusive_target_contrast.csv - Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).
stagewise_function\per_regulon_target_overlap.csv - Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).
stagewise_function\target_coverage.csv - Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).

Definitions

Adult changes are relapse minus diagnosis, centered over 201 background regulons and
oriented to the pediatric effect. Positive values denote agreement with that direction.
Null A flips patient-level signs; null B uses size-matched random targets; null C draws
55-regulon modules from the 201-regulon background with their pediatric orientations.
The empirical p-value denominator is 2,001 for 2,000 draws. Per-regulon paired p values
are unadjusted. Functional-enrichment q values use the families specified in Methods.
Pooled enrichment separates positive and negative TF-gene correlation signs. These signs
are distinct from pediatric orientations and do not establish activation or repression.
Marker-excluded module scores retain the complete maturity axes. Greater axis values
indicate primitive marker expression. The earlier weakened-axis sensitivity is separate.
Stage-exclusive enrichment tables include nonsignificant tests and the common-background
contrast, so differences between separate lists are not treated as a pathway-switching test.
