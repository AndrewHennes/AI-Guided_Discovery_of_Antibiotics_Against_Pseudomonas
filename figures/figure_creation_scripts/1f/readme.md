# Figure 1f generating code

`plot_growth_inhibition_distribution.py` contains the data processing and plotting implementation used for the final version 025 panel. Run this file directly with the Python packages in the methods root's `requirements.txt`. No command-line arguments are needed, and all paths resolve absolutely from the script location.

The script reads the bundled source tables from `shared/figure1_current/inputs` and writes its plot, supporting tables, and validation to `1f/generated_from_data`. This separate output directory preserves the supplied final artwork. The final Illustrator document also contains its native typography and panel placement, which are retained by the existing artwork-export pipeline.

The script includes the complete fixed-bandwidth Gaussian density calculation and shared logarithmic width transform. Its calculated curves are checked against the verified final-panel tables. The published panel f artwork is preserved.

See `figure_1_panel_source_code.md` at the methods root for the panel-to-code map and provenance.

## Existing artwork export

# Figure 1f

Checks the continuous measurements and recalculated logarithmic-density curves against the verified tables, then exports the unchanged native panel f. Screen n = 7,994 is retained.

Run `growth_inhibition_log_density_violins_f.py` to regenerate this panel and synchronize its SVG in `figures/figure_svg_files`. Paths resolve absolutely from the script location, so the working directory does not matter. The script takes no command-line arguments. See `README_Figure_1.md` at the methods root for provenance and validation.
