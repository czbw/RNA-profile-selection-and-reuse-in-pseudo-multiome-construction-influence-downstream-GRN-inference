Supplementary Data 1. Published link and metacell comparisons

CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.
The table index identifies source files and selection rules. Source files are preserved under source data/.
source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.
No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.

Tables:
membership.csv - TF and signed eRegulon membership (Figure 1C).
shared_regulon_agreement.csv - Activity and target agreement for the same shared eRegulons (Figure 1D).
shared_tf_targets.csv - Target overlap after pooling signed eRegulons within each TF.
network_counts.csv - Network sizes underlying Figure 1C.

Definitions

The two signs in an eRegulon identity refer to TF-gene and region-gene correlations.
Activity agreement is Pearson correlation across cell-type mean AUCell profiles.
Target-gene and target-region Jaccard indices refer to the same shared eRegulons.
