# Figure 1g generating code

`plot_all_hit_properties.py` contains the data processing and plotting implementation used for the final version 025 panel. Run this file directly with the Python packages in the methods root's `requirements.txt`. No command-line arguments are needed, and all paths resolve absolutely from the script location.

The script reads the bundled source tables from `shared/figure1_current/inputs` and writes its plot, supporting tables, and validation to `1g/generated_from_data`. This separate output directory preserves the supplied final artwork. The final Illustrator document also contains its native typography and panel placement, which are retained by the existing artwork-export pipeline.

`amine_descriptor.py` in this folder contains the exact primary, secondary, and tertiary amine definition and 33 functional-group checks. All-hit groups contain WT 8, Δ40 5, Δ3 39, and Δ6 124 compounds, including shared hits in each relevant group.

See `figure_1_panel_source_code.md` at the methods root for the panel-to-code map and provenance.

## Existing artwork export

# Figure 1g

Recalculates TPSA, primary + secondary + tertiary amine count, and LogP for all active compounds in WT, Δ40, Δ3, and Δ6. Counts are 8, 5, 39, and 124. Shared hits are retained in each applicable group.

Run `tpsa_amine_count_log_p_all_strain_hits_g.py` to regenerate this panel and synchronize its SVG in `figures/figure_svg_files`. Paths resolve absolutely from the script location, so the working directory does not matter. The script takes no command-line arguments. See `README_Figure_1.md` at the methods root for provenance and validation.
