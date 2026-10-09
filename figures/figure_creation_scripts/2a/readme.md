# Figure 2a

Run `deepflux_gi60_gi80_delta5_hits_vs_inactive_roc_a.py` to reproduce the ROC plot. The filtered CO-ADD dataset supplies one observation per `COADD_Compounds_ID`. Δ5 hits are positive; compounds inactive against both WT and Δ5 are negative. Both the GI60 and GI80 models are evaluated against the 80% growth-inhibition hit labels. Missing measurements are excluded.

The cohort contains 76 positive and 32,228 negative compound IDs. One positive ID is active against both strains. Labels, standardized structures, and ensemble scores agree among the eligible records belonging to each ID. The first eligible source record is retained. Distinct compound IDs remain separate even when they share a molecular structure. `ROC_compound_ID_record_counts.csv` records how many eligible input rows belong to each ID and identifies the representative source row.

The analysis uses the original GI60 and GI80 ensemble means, giving AUCs of 0.7691685230 and 0.7502314921, respectively. False positive rate is horizontal and true positive rate is vertical. All empirical ROC coordinates are saved; collinear vertices are omitted from the rendered paths.

Run `score_coadd_with_deepflux_gi60_gi80_a.py` to recreate the full scored dataset from the included checkpoints and MiniMol vectors. `ROC_Eligible` marks records with eligible measured labels. `ROC_Included` marks the single selected observation per compound ID. `ROC_Eligible_Records_For_ID` gives the eligible record count for each ID. The original probabilities and hit labels are retained. The selected cohort is saved as `COADD_Delta5_hits_and_inactive_ROC_cohort.csv`.

Dataset paths are blank in `../figure2_config.py`. Saved-mode reproduction reads all 34,661 eligible input records from `saved_inputs/COADD_roc_input_scores.csv.gz`, validates within-ID agreement, and selects 32,304 unique compound IDs. This removes 2,357 extra eligible records. The scoring script refreshes the same saved input and cohort. The detailed analysis definition is in `ROC_definition.json`.
