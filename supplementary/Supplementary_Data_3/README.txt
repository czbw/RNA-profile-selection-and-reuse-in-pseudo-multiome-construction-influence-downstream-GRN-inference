Supplementary Data 3. Matched-reuse assignments and activity

CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.
The table index identifies source files and selection rules. Source files are preserved under source data/.
source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.
No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.

Tables:
matched_reuse_summary.csv - Assignment counts and activity variance for all six matched-reuse arms (Figure 3A,D).
multiplicity_checks.csv - Per-stratum equality checks for the two reuse-matched comparisons.
stratum_census.csv - The 190 donor-cell-type strata used in matched-reuse comparisons.
donor_activity_summary.csv - Donor-stratum mean activity agreement (Figure 3C), with eligibility sensitivity results.
donor_activity_by_regulon.csv - Per-cell-type-regulon donor-stratum correlations for the primary eligibility rule.
capacity_sweep.csv - Capacity sweep used to choose the intermediate reuse level.

Definitions

Reuse ratio is the number of ATAC assignments divided by distinct assigned RNA profiles.
All six arms contain 69,249 assignments. The three levels use 69,249, 49,262 and 34,488
distinct RNA profiles. A stratum is a donor-label and cell-type combination. At reuse
1.406 and 2.008, all 190 strata have equal inferred and random multiplicity distributions
(KS distance zero); one-to-one arms use each eligible RNA profile once by construction.
Variance ratios have the native-pairing AUCell variance as denominator, not RNA count variance.
rho_strata is the median Spearman correlation over eligible cell-type-regulon combinations.
The primary rule requires 20 cells per stratum and five donor strata per cell-type.
Suffixes c10_d5, c50_d5, c20_d3 and c20_d7 give alternative eligibility thresholds.
