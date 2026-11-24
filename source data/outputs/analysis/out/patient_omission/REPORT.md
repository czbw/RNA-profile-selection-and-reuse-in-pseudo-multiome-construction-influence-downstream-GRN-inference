# Patient-omission results - 29 September 2026

The overall expression-score signal survives every single-patient omission, while exact nomination membership is substantially more sensitive. Keep the current 55 as the recorded retrospective set; do not silently replace it with a selected subset.

## Nomination stability

- Full-data reconstruction reproduced all three axes and the exact 55-member set. Numerical comparisons of E2 effects, compartment effects, cell-type counts, stored orientations, external patient deltas and pediatric decomposition passed.
- Across 19 omitted pediatric patients, newly selected sets contain 43–60 regulons; 39–53 of the original 55 remain.
- 26/55 survive every omission; 29/55 survive at least 90% (at least 18/19) of omissions. The least stable original member survives only 2/19.
- Jaccard with the original set ranges from 0.606 to 0.898.
- Recomputing the at-least-five-patient compartment eligibility sets versus holding those sets fixed gives identical membership results. Omitting the healthy-only patient also leaves membership unchanged.

The 26 members present in every omission are a descriptive stability subset. They are not the historical 26-regulon module or the corrected LMM-significant 26, and are not proposed as a replacement module.

### Least stable original members

| Regulon | Selected folds | Ax1 passes | Ax2 passes | Ax3 passes |
|---|---:|---:|---:|---:|
| EP300_direct_-/- | 2/19 | 19 | 10 | 11 |
| EP300_direct_-/+ | 7/19 | 19 | 15 | 11 |
| GABPB1_direct_-/+ | 7/19 | 19 | 13 | 12 |
| KLF10_direct_+/+ | 10/19 | 18 | 17 | 12 |
| BACH1_direct_+/+ | 11/19 | 18 | 12 | 18 |
| HMGA1_direct_+/+ | 11/19 | 19 | 11 | 19 |
| KLF10_direct_+/- | 11/19 | 14 | 18 | 11 |
| MAFB_direct_+/- | 11/19 | 19 | 11 | 19 |
| KLF3_direct_+/+ | 12/19 | 19 | 12 | 19 |
| KLF4_direct_+/- | 12/19 | 19 | 12 | 19 |
| MAFF_direct_+/- | 12/19 | 12 | 19 | 12 |
| ERG_direct_-/- | 13/19 | 19 | 13 | 19 |

Axis 3 was historically evaluated only for axis-1 survivors; these pass counts are not independent failure categories. Complete decisions for all 201 candidates and every fold are retained.

## Fixed-module patient influence

| Cohort / metric | Targets | Full mean | Omission mean range | Largest omission Wilcoxon p |
|---|---|---:|---:|---:|
| external / module_delta | Diagnosis | 0.034375 | 0.032146–0.036610 | 1.49e-07 |
| external / module_delta | Relapse | 0.025516 | 0.023846–0.027309 | 6.29e-07 |
| pediatric / composition | Diagnosis | 0.002140 | -0.000070–0.004185 | 0.64 |
| pediatric / composition | Relapse | 0.001717 | -0.000072–0.003435 | 0.64 |
| pediatric / equal_celltype_weight | Diagnosis | 0.009754 | 0.008739–0.010974 | 0.00336 |
| pediatric / equal_celltype_weight | Relapse | 0.008069 | 0.006999–0.009030 | 0.0056 |
| pediatric / total | Diagnosis | 0.012465 | 0.010121–0.014852 | 0.0304 |
| pediatric / total | Relapse | 0.010293 | 0.009004–0.012527 | 0.0139 |
| pediatric / within_type | Diagnosis | 0.010325 | 0.008806–0.011712 | 0.00233 |
| pediatric / within_type | Relapse | 0.008577 | 0.007542–0.009637 | 0.00281 |

The overall, equal-cell-type-weighted and within-type module effects retain their positive sign after every omission. The smaller composition component can change sign and its sensitivity p-values are not consistently below 0.05. The decomposition uses six admitted types and diagnosis weights; it is not an edge-rewiring measurement or whole-marrow decomposition.

## Re-selection and orientation sensitivity

- Diagnosis targets: external mean after pediatric fold-wise re-selection/orientation ranges 0.026504–0.036824; largest Wilcoxon p=1.43e-07. The omitted-patient score has mean 0.008523, 78.9% positive and descriptive Wilcoxon p=0.00945.
- Relapse targets: external mean after pediatric fold-wise re-selection/orientation ranges 0.019727–0.027798; largest Wilcoxon p=2.57e-06. The omitted-patient score has mean 0.006421, 78.9% positive and descriptive Wilcoxon p=0.00823.

These omitted-patient scores avoid using that patient in the downstream membership/sign calculation. They are not independent validation: network inference, target definitions and score panels still use the full discovery cohort. Fold scores also share training patients, so their descriptive p-values should not be read as a confirmatory cross-validation test.

## Implementation and interpretation

- `f6_10_patient_omission.py` uses existing pseudobulk and external AUCell matrices. Pediatric centering uses the complete stage-specific panel (332/352 columns), whereas external centering uses its recorded 201-regulon background. Subsetting before centering would change the analysis.
- Six cell types feed directional agreement and the cell-type rule; seven malignant and eight healthy types feed the compartment comparison. Both stages and compartments are removed together for a patient.
- All target sets and networks stay fixed; no CUBIC execution is needed. This does not measure uncertainty in network inference or establish generalization to another cohort.
- No original results, module definitions, thresholds or main figures were replaced. Full inputs and software hashes are in `run_manifest.json`; full-cohort parity is in `baseline_parity.json`.
- The frozen external orientation is the sign of E2 effects: 30 positive and 25 negative module members under both vocabularies. The funnel’s displayed 31-up/24-down split uses a different direction field derived from the compartment difference. These must not be conflated.

Suggested Results wording: “Patient-omission analysis supported the aggregate expression-score association but revealed sensitivity of exact membership: 26 of 55 nominated regulons were retained in all 19 omissions, and each fold retained 39–53 members. External module-score changes remained positive under both fixed-set patient deletion and discovery-fold re-selection. These analyses were conditional on the original inferred networks and target definitions.”
