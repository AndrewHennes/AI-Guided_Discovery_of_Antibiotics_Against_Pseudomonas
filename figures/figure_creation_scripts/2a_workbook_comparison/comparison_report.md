# CO-ADD ROC comparison and repeated-entry audit

One Δ5-hit structure contributes 64 CSV records split equally between CO-ADD:0307276 and CO-ADD:0302755. Its GI60 and GI80 probabilities are approximately 0.001845 and 0.001797. This score pair occurs twice in the workbook. The analyses below compare record, compound-ID, and structure counting units.

## Earlier workbook and unfiltered CSV

| Analysis at an 80% hit cutoff | Positive records | Negative records | GI60 AUC | GI80 AUC |
|---|---:|---:|---:|---:|
| Earlier workbook, both models evaluated against its 80% labels | 92 | 34,106 | 0.776213 | 0.740836 |
| Filtered CSV, every eligible row | 143 | 34,518 | 0.594723 | 0.594668 |
| Unfiltered CSV, every eligible row | 154 | 34,630 | 0.605315 | 0.592573 |
| Unfiltered CSV, reducing only the 64-row hit to its two compound IDs | 92 | 34,630 | 0.776840 | 0.740954 |
| Unfiltered CSV, one row per CO-ADD compound ID throughout | 87 | 32,336 | 0.765997 | 0.726989 |
| Unfiltered CSV, one row per standardized structure throughout | 86 | 32,287 | 0.770801 | 0.731100 |

The single-hit diagnostic reduces that structure to its two compound IDs. The global analyses deduplicate the full cohort.

The workbook originally evaluates GI60 against `Binarized_at_60_PGI`, giving 155 positives, 34,043 negatives, and AUC 0.768410540. It evaluates GI80 against `Binarized_at_80_PGI`, giving 92 positives, 34,106 negatives, and AUC 0.740836433. Both original areas exactly reproduce the associated Prism project. Evaluating both models at an 80% cutoff gives the workbook values in the table above.

The stored score pairs for 32,835 of the workbook's 34,198 records closely match scores in the newly scored unfiltered CSV, using relative tolerance 0.00001 and absolute tolerance 0.000000006 in each score. Within these numerical matches, recomputed probabilities give effectively identical AUCs to the workbook probabilities when holding workbook labels fixed.

## GI60 evaluated at 60% inhibition

CSV positives have Δ5 inhibition ≥60%. Negatives have inhibition <60% in both Δ5 and WT. Missing measurements are excluded. All Δ5-negative records with measured Δ5 values in these CSV cohorts also satisfy the WT-negative condition at this cutoff. Workbook results use its supplied 60% binary label.

| Dataset | Counting unit | Positive records | Negative records | GI60 AUC |
|---|---|---:|---:|---:|
| Earlier workbook | Workbook record | 155 | 34,043 | 0.768411 |
| Filtered CSV | Every eligible row | 205 | 34,456 | 0.644796 |
| Filtered CSV | One per CO-ADD ID | 136 | 32,168 | 0.763384 |
| Filtered CSV | One per standardized structure | 135 | 32,119 | 0.766410 |
| Unfiltered CSV | Every eligible row | 217 | 34,567 | 0.649749 |
| Unfiltered CSV | One per CO-ADD ID | 148 | 32,275 | 0.761105 |
| Unfiltered CSV | One per standardized structure | 147 | 32,226 | 0.763870 |
| Unfiltered CSV | Only the 64-row hit reduced to its two IDs | 155 | 34,567 | 0.769134 |

The Δ5 binary labels agree within each CO-ADD ID and standardized structure at both the 60% and 80% cutoffs. The reproduction script verifies this before selecting one row from each group.

## Repeated entries in the other datasets

The following results refer to the filtered datasets used to recreate the figures. A repeated structure means the same canonical isomeric SMILES occurs on multiple rows. Extra rows are records after the first for each structure. Exact duplicate rows match every original column.

| Dataset | Rows | Repeated structures | Extra rows for repeated structures | Exact extra duplicate rows |
|---|---:|---:|---:|---:|
| Broad 8K | 7,997 | 12 | 12 | 0 |
| Broad 39K | 39,300 | 543 | 1,073 | 607 |
| ECBD | 100,763 | 0 | 0 | 0 |
| CO-ADD | 185,502 | 2,476 | 86,404 | 0 |

All repeated structures in the filtered 8K dataset are inactive for all four strains. Its 124 Δ6 hits, including 8 WT hits, are unique structures, so the Figure 2d Gaussian input cohort is unaffected. Its Figure 1 strain hit counts are also unaffected.

The unfiltered 8K dataset has 8 repeated structures and 8 extra rows, with one repeated Δ6/Δ3 hit. Unfiltered ECBD has no repeated structures. Unfiltered CO-ADD has 2,495 repeated valid structures and 89,002 extra rows. Source rows with invalid or missing structures are reported separately and excluded from these structure counts. The CO-ADD rows differ in other assay or identifier columns even when they repeat the same Δ5 observation, so their full-row exact-duplicate count is zero.

The complete seven-dataset audit, repeated-row extracts, per-strain label agreement, and reproduction script are in `outputs/dataset_repeat_audit`. Dataset inputs were not modified.
