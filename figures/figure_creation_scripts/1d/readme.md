# Figure 1d generating code

`plot_dataset_hit_counts.py` contains the data processing and plotting implementation used for the final version 025 panel. Run this file directly with the Python packages in the methods root's `requirements.txt`. No command-line arguments are needed, and all paths resolve absolutely from the script location.

The script reads the bundled source tables from `shared/figure1_current/inputs` and writes its plot, supporting tables, and validation to `1d/generated_from_data`. This separate output directory preserves the supplied final artwork. The final Illustrator document also contains its native typography and panel placement, which are retained by the existing artwork-export pipeline.

See `figure_1_panel_source_code.md` at the methods root for the panel-to-code map and provenance.

## Existing artwork export

# Figure 1d

Recalculates unique active structures for the supplied 8K, 39K, ECBD, and CO-ADD inputs and rebuilds the count plot.

Run `unique_hit_counts_across_datasets_d.py` to regenerate this panel and synchronize its SVG in `figures/figure_svg_files`. Paths resolve absolutely from the script location, so the working directory does not matter. The script takes no command-line arguments. See `README_Figure_1.md` at the methods root for provenance and validation.

## Generation code

The underlying plotting and analysis code is collected in [source_code](source_code/readme.md), with its input references and a reproduction entry script. The current panel export script above retains the final version 025 arrangement.
