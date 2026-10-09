# Figure 2d

Run `gaussian_property_distributions_delta6_only_and_wild_type_hits_d.py` to calculate the properties, fit Gaussian distributions, and export the plot. The cohort contains 116 Δ6-only hits and eight WT/Δ6 hits from the 8K dataset.

The output table contains the source identifiers, SMILES, binary hit calls, category, IPC, log10IPC, logP, and number of aromatic rings. Separate CSV files contain the fitted parameters, 10th/25th/50th/75th/90th percentiles, and two-sided Welch test results.

Dataset paths are blank in `../figure2_config.py`. The default saved mode uses the included 8K input table. See `../README_Figure_2.md` for the property definitions and plotting settings.

Run `compare_recreated_and_supplied_gaussian_properties_d.py` to compare these fits with the supplied efflux-rule Gaussian parameters. The comparison, individual plots, parameters, and percentile differences are written to `gaussian_comparison`. The supplied parameter columns map Hits to WT/Δ6 active and Nonhits to Δ6-only active; IPC is treated as log10IPC. Both panels use the same axes, fonts, colors, box widths, and percentile markers. Their common vertical range is −0.5 to 8, including the supplied WT aromatic-ring 10th percentile. Comparison panels display the distributions without significance annotations. Class priors are recorded separately and do not affect within-class Gaussian percentiles.
