Supplementary Data 6. Pediatric decomposition and sensitivity analyses

CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.
The table index identifies source files and selection rules. Source files are preserved under source data/.
source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.
No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.

Tables:
module_decomposition.csv - Pediatric module decomposition; oriented centered scores are primary (Figure 5C).
patient_decomposition.csv - Patient-level total, within-cell-type and composition components.
regulon_decomposition.csv - Per-regulon pediatric decomposition and statistical comparisons.
clp_sensitivity_summary.csv - Seven-cell-type sensitivity including CLP.
clp_candidate_criteria.csv - CLP sensitivity results.
clp_regulon_decomposition.csv - CLP sensitivity results.
clp_adult_patient_changes.csv - CLP sensitivity results.
clp_adult_scores.csv - Adult score and null summaries for the CLP sensitivity.
clp_fixed_module_adult_scores.csv - Adult score and null summaries for the CLP sensitivity.
patient_omission\all_candidate_fold_decisions.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\baseline_rebuilt_nomination.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\external_sample_pairs.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\fixed_module_patient_deltas.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\fixed_module_patient_influence.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\nomination_by_omitted_patient.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\regulon_selection_frequency.csv - Patient-omission results (Supplementary Figures S1,S2).
patient_omission\reselected_module_scores.csv - Patient-omission results (Supplementary Figures S1,S2).

Definitions

oriented_offsetcentred is the primary score scale. offsetcentred is the secondary
unoriented centered scale. raw is an additional uncentered, unoriented scale.
Total change equals within-cell-type change at diagnosis composition plus the residual
composition component. The composition share can be unstable when total change is near zero.
Source columns containing rewiring denote within-cell-type activity change, not network edges.
Patient-omission analyses condition on the inferred networks and do not refit them.
