# Figure 2

Run `reproduce_figure_2.py` to calculate panels a and d, export panels b, c, e, and f, and assemble Figure 2. Each panel has its own folder, script, PDF, SVG, and PNG. The two unlabeled lower panels are in `2_unlabeled_growth_inhibition` and `2_unlabeled_compound_examples`. The complete figure is written to `full_figure_2`.

The default `DATA_MODE = "saved"` uses the included input tables. Dataset paths are blank in `figure2_config.py`. To use the external datasets, set `DATA_MODE = "datasets"` and enter absolute CSV paths. All output paths are resolved from the scripts' absolute locations.

The vector layout supplies the original fonts, axes, labels, marker geometry, and artwork. Panels a and d replace the plotted values with calculations from the supplied data. The remaining panels retain their plotted geometry and use the Figure 1 strain colors. Panel e imports the supplied scatter SVG components, with open black prediction circles matching the green WT markers in size. The six embedded background images and individual vector markers are retained. The model-series colors in panel a and the compound-identity symbols remain distinct from the strain palette.

The PDF retains embedded fonts and text. The SVG uses vector text outlines. Run `save_native_illustrator_figure_2.jsx` from Illustrator's File → Scripts → Other Script menu to import the assembled PDF and save a fresh native Illustrator file with PDF compatibility enabled.

## Panel a

The receiver operating characteristic (ROC) analysis compares Δ5 hits with compounds inactive against both WT and Δ5 in the filtered CO-ADD dataset. A positive record has `PA_MexB_Binarized_at_80_Percent_Growth_Inhibition = 1`, including hits that are also active against WT. A negative record has both `PA_MexB_Binarized_at_80_Percent_Growth_Inhibition = 0` and `PA_Binarized_at_80_Percent_Growth_Inhibition = 0`. Binary hits use the 80% growth-inhibition threshold. Missing measurements are not assigned to the inactive class.

The cohort contains 76 positive and 32,228 negative CO-ADD compound IDs. Of the positives, 75 IDs are Δ5-only hits and one is active against both strains. The positive and negative groups contain 75 and 32,179 distinct standardized SMILES, respectively. Each `COADD_Compounds_ID` contributes once. Eligible source records are checked for identical hit labels, standardized structures, and ensemble means within each ID, and the first eligible record is retained. The 34,661 eligible source records yield 32,304 compound IDs, removing 2,357 additional records. All selected IDs have scores for both models.

GI60 and GI80 are both evaluated against the 80% inhibition cutoff. Scores are the arithmetic means of sigmoid probabilities from the models' ten ensemble members. Larger original ensemble scores rank observations toward the Δ5-hit positive class. The areas under the ROC curve are 0.7691685230 for GI60 and 0.7502314921 for GI80. False positive rate is horizontal and true positive rate is vertical. All empirical score thresholds are saved in the coordinate table; collinear vertices are omitted from rendered paths. Legends display AUC to three decimal places.

`2a/score_coadd_with_deepflux_gi60_gi80_a.py` recreates the scored CO-ADD table using the included checkpoints and 512-dimensional MiniMol representations. Inference runs on the CPU in evaluation mode. Each checkpoint is loaded with `weights_only=True`. Per-member scores, ensemble means, population standard deviations, source row indices, cohort labels, and representation indices are included. The inference implementation was checked against the supplied model class for every checkpoint, with identical logits for 256 input vectors per checkpoint.

The representation archive, `2a/saved_inputs/COADD_minimol_features.npz`, stores `smiles` and `features` arrays in the same order. `Minimol_Feature_Index` is the zero-based row in that archive. The source table contains 185,502 records, of which 185,494 have scores.

The complete scored table retains every input record. `ROC_Eligible` identifies records satisfying the class definitions, and `ROC_Included` identifies the one selected observation per ID. `ROC_Eligible_Records_For_ID` records each ID's eligible source-record count. `ROC_compound_ID_record_counts.csv` lists the record count and representative source row for every included ID.

`COADD_Delta5_hits_and_inactive_ROC_cohort.csv`, `ROC_metrics.csv`, and `ROC_curve_coordinates.csv` reproduce the plotted comparison. The label is `ROC_Delta5_Hit_Positive_Label`; ranking scores are `ROC_GI60_Delta5_Hit_Score` and `ROC_GI80_Delta5_Hit_Score`, equal to their original ensemble means. `saved_inputs/COADD_roc_input_scores.csv.gz` retains all eligible source records so the per-ID selection is reproduced when the script runs. External dataset paths remain blank. Checkpoint checksums and inference settings are in `scoring_report.json`.

## Panel d

The 8K cohort contains 124 Δ6 hits. The WT and Δ6 active category contains the eight records with both binary hit columns equal to 1. The Δ6-only category contains the 116 records with the Δ6 binary column equal to 1 and the WT binary column equal to 0.

Properties are calculated from `Standardized_SMILES` with RDKit. `IPC` is `Descriptors.Ipc` with its default arguments. `log10IPC` is its base-10 logarithm. `logP` is `Crippen.MolLogP`. `Number_of_Aromatic_Rings` is `CalcNumAromaticRings`.

A Gaussian distribution is fit separately to each property and category by maximum likelihood. The fitted standard deviation uses the population denominator. Boxes span the fitted 25th through 75th percentiles, solid black lines mark the fitted 50th percentile, and capped whiskers mark the fitted 10th and 90th percentiles. Percentiles are calculated from the fitted normal distribution without random sampling.

The three asterisks denote a two-sided Welch t-test with p < 0.05 on the observed property values. Exact statistics and degrees of freedom are saved in `Property_comparison_statistics.csv`. The labeled molecule table and Gaussian parameters are saved beside the plot.

## Scripts

| Folder | Script |
| --- | --- |
| `2a` | `deepflux_gi60_gi80_delta5_hits_vs_inactive_roc_a.py` |
| `2a` | `score_coadd_with_deepflux_gi60_gi80_a.py` |
| `2b` | `ranked_compound_delta6_and_wild_type_hit_retrieval_b.py` |
| `2c` | `compound_ic50_in_wild_type_and_delta6_strains_c.py` |
| `2d` | `gaussian_property_distributions_delta6_only_and_wild_type_hits_d.py` |
| `2e` | `deepflux_efflux_score_compound_predictions_and_wild_type_hits_e.py` |
| `2f` | `molecular_property_distributions_by_growth_inhibition_f.py` |
| `2_unlabeled_growth_inhibition` | `growth_inhibition_cumulative_distributions_figure_2_unlabeled.py` |
| `2_unlabeled_compound_examples` | `compound_structures_and_dose_responses_figure_2_unlabeled.py` |
| Methods root | `reproduce_figure_2.py` |

Use the packages pinned in `requirements.txt`; Figure 2 inference additionally uses the pinned `torch` package. The included MiniMol representations are the inference inputs, so regenerating molecular embeddings is not required to reproduce the scores or plots.

Additional CO-ADD analyses are saved in `2a_unfiltered` and `2a_workbook_comparison`. The first scores the unfiltered CO-ADD dataset and plots its ROC curves. The second compares the earlier workbook, 60% and 80% hit cutoffs, and record-counting units. Each folder has its own reproduction scripts and README; they export separate comparison figures. `2a_workbook_comparison/comparison_report.md` records the repeated-entry audit findings and numerical results.
