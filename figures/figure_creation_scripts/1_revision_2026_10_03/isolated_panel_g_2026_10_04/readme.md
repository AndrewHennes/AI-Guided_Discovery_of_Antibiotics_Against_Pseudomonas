# Standalone property panels requested on 4 October 2026

`reproduce_isolated_properties.py` creates exactly three SVG files, a combined PNG preview, and the supporting tables. All input and output paths are absolute. There are no command-line arguments. The script runs with the established Chemprop Python environment at `/opt/anaconda3/envs/Chemprop_Environment/bin/python`.

The script reads the current supplied 8K CSV directly, recomputes log10 IPC, logP, and aromatic ring count for the union of hit structures, and creates four groups using the strain-specific binary activity labels. Counts are WT 8, Δ40 5, Δ3 39, and Δ6 124. Groups overlap. All 116 Δ6 hits inactive in WT remain included within the full Δ6 group. No inactive comparison group is drawn.

All compounds are shown as points in every group. The point layer includes outliers, so box-plot outlier markers are suppressed to avoid duplication. Each SVG is an editable vector panel with no embedded raster image. The current main Illustrator figure and its earlier reproduction scripts are not edited by this workflow.

`validate_isolated_properties.py` verifies the exported descriptor values, strain membership, and actual point counts inside each SVG. Definitions and caption text are in the accompanying `METHODS_AND_CAPTION.md`. The table `hit_group_overlap.csv` records shared compounds, and `plotted_compounds.csv` lists every plotted observation.
