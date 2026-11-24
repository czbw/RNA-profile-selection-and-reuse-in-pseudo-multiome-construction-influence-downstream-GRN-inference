Supplementary Data 4. Network agreement and pAML nomination

CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.
The table index identifies source files and selection rules. Source files are preserved under source data/.
source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.
No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.

Tables:
matched_reuse_networks.csv - Absolute agreement with the native network (Figure 3B and Supplementary Figure S3).
paml_reuse.csv - pAML RNA reuse using sample-plus-barcode identities (Figure 4).
paml_network_counts.csv - Published link and corrected network sizes (Figure 4).
corrected_nomination.csv - All shared candidates and the three further cumulative nomination criteria (Figure 4).
corrected_paired_tests.csv - pAML patient-cell-type units and longitudinal statistical results.
corrected_mixed_models.csv - pAML patient-cell-type units and longitudinal statistical results.
corrected_paired_units.csv - pAML patient-cell-type units and longitudinal statistical results.
corrected_calibration.csv - Mixed-model calibration diagnostics; numeric results only.
published_link_nomination.csv - All shared candidates and the three further cumulative nomination criteria (Figure 4).
published_link_paired_tests.csv - pAML patient-cell-type units and longitudinal statistical results.
published_link_mixed_models.csv - pAML patient-cell-type units and longitudinal statistical results.
published_link_paired_units.csv - pAML patient-cell-type units and longitudinal statistical results.
published_link_calibration.csv - Mixed-model calibration diagnostics; numeric results only.
nomination_counts.csv - Cumulative candidate counts (Figure 4).
paml_stratum_census.csv - pAML eligible RNA and ATAC counts by patient, stage and cell-type.
paml_stratum_assignment.csv - pAML eligible RNA and ATAC counts by patient, stage and cell-type.

Definitions

Shared signed eRegulon identity is the first nomination criterion. The three further
criteria are consistent longitudinal direction, malignant-versus-healthy selectivity and
cell-type agreement. nomination_counts applies these criteria cumulatively. Retired
TF-family and ChIP criteria have not been carried into the curated candidate tables.
pAML reuse counts distinguish RNA profiles by sample and barcode, not barcode alone.
The older source field named distinct_donors refers to RNA profiles, not clinical donors.
