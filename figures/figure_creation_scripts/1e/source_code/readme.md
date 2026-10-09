# Figure 1e generation code

The binary hit table and the code that calculates the strain intersections and updates the Venn counts.

`shared/analysis.py` and `shared/rendering.py` are copies of the code used for the source figure. The input layout and saved plot tables are in the methods root's `shared` directory. `source_inputs.json` identifies them by checksum. `original_config.py` records the original settings; `config.py` retains those settings and sets absolute locations for this installation.

Run `reproduce_source_1e.py` to reproduce the source panel and its analysis tables. Results are written to `/Users/asselism/Documents/Codex/2026-09-21/please-create-an-html-file-that/outputs/final_figure1_source_reproduction/1e`. External dataset paths remain blank. The scripts have no command-line arguments and can run from any working directory.

The current Figure 1 positions and labels are exported by `../four_strain_antibacterial_hit_overlap_venn_e.py`. That existing entry script uses the version 025 layout. The source reproduction keeps the earlier source-panel layout and label e, from which the current panel was assembled.
