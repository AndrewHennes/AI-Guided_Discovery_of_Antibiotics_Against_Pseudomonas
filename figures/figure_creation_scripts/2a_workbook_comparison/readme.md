# CO-ADD workbook and growth-cutoff comparison, Figure 2a

Run `coadd_workbook_and_growth_cutoff_comparison_2a.py` to regenerate all comparison plots, ROC metrics, and plotted coordinates. It reads the saved workbook values here and the scored inputs in `../2a` and `../2a_unfiltered`. The optional external input paths are blank. All paths resolve absolutely from the script location.

The CSV analysis uses Δ5 inhibition at or above the selected cutoff as positive and inhibition below the cutoff in both Δ5 and WT as negative. The workbook analysis uses its supplied binary columns. The original ensemble probabilities are used throughout. Curves use false positive rate on the horizontal axis and true positive rate on the vertical axis. AUCs are checked independently using ranks and curve integration.

The saved workbook columns named `GI60_Hits` and `GI60_NonHits` exactly match the GI60 scores split by `Binarized_at_60_PGI`. Likewise, `GI80_Hits` and `GI80_NonHits` match GI80 scores split by `Binarized_at_80_PGI`.

`GI60_60_percent_cutoff_comparison_2a.svg` compares the earlier workbook, filtered CSV records, unfiltered CSV records, and one observation per unfiltered CO-ADD ID. `workbook_and_unfiltered_80_percent_comparison_2a.svg` compares the workbook at a common 80% cutoff with the unfiltered CSV before and after adjusting the single repeated hit. Full numerical results are in `ROC_comparison_metrics.csv`.

See `comparison_report.md` for the numerical comparison and the separate dataset repeated-entry audit in the outputs directory.
