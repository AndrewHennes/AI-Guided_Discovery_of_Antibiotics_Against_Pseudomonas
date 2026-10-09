# Figure methods

Figure 1 uses the current version 025 layout, with panels **a–i**. Run `reproduce_figure_1.py` to recompute the current counts and properties, verify the inhibition distributions, and refresh every Figure 1 SVG in `../figure_svg_files`. See [README_Figure_1.md](README_Figure_1.md) for the current panel map, inputs, exports, and validation. [METHODS_Figure_1.md](METHODS_Figure_1.md) supplies methods and caption text.

Figure 2 instructions remain in [README_Figure_2.md](README_Figure_2.md). Its scripts, data, and SVG exports are unchanged by the Figure 1 update.

The standalone dataset comparison and offline molecule explorer remain in [dataset_hit_counts](dataset_hit_counts/README.md). Additional analyses remain in [undesigned](undesigned/README.md).

The previous Figure 1 a–g pipeline and SVGs are preserved in `archive/figure_1_before_v025_sync_2026_10_05.zip`. The folder `1_revision_2026_10_03` preserves the development history and native Illustrator editing scripts. Use the current root entry script for publication exports. Root `config.py` and `shared/rendering.py` support historical artwork and Figure 2 utilities; current Figure 1 settings are in `figure1_current_config.py`.

The older shared validation reports describe the export revisions named in those reports. Current Figure 1 validation is recorded in `figure_1_current_validation.json` and `current_figure_1_data/analysis_validation.json`.

## Environment

Use Python 3.11 or newer with `requirements.txt`. Scripts resolve inputs and outputs from their own absolute paths, so they can run from any working directory. No command-line arguments are needed.
