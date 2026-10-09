# Unfiltered CO-ADD ROC analysis, Figure 2a

Run `deepflux_unfiltered_coadd_delta5_hits_vs_inactive_roc_2a.py` to regenerate the SVG, PDF, PNG, metrics, and ROC coordinates from the included scored records. The script uses absolute paths resolved from its own location. `SCORED_DATASET_PATH` is blank and can be set to an absolute path for a different scored input.

The positive class is Δ5 growth inhibition ≥80%. The negative class is growth inhibition <80% in both Δ5 and WT. Missing measurements are excluded. Every eligible input record contributes once. Both ensembles use their original mean probabilities.

`score_unfiltered_coadd_with_deepflux_2a.py` regenerates the full scored dataset from the saved input, the feature caches in this folder and `../2a`, and the 20 checkpoints in `../2a/checkpoints`. Its optional external dataset path is blank. CPU inference uses 8 threads. Install the dependencies listed in the parent methods directory before running the scripts.

The saved input contains 188,740 records, with scores for 188,730. At the 80% cutoff, 154 positives and 34,630 negatives yield AUCs of 0.6053148122 (GI60) and 0.5925726886 (GI80).

These exports are an additional analysis. The main Figure 2 exports remain separate.
