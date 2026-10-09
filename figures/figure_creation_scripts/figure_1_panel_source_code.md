# Generating code for the final Figure 1 panels

The supplied `Figure1_v025_all_hits_TPSA_amines_LogP (3).ai` is byte-identical to the verified version 025 native figure. The full plotting code for the panels created in this revision is now included in the corresponding panel directories.

| Final panel | Plot | Generating code |
| --- | --- | --- |
| d | Dataset hit counts | [1d/plot_dataset_hit_counts.py](1d/plot_dataset_hit_counts.py) |
| f | Growth-inhibition log-density violins | [1f/plot_growth_inhibition_distribution.py](1f/plot_growth_inhibition_distribution.py) |
| g | All-hit TPSA, amine count, and LogP boxplots | [1g/plot_all_hit_properties.py](1g/plot_all_hit_properties.py) |
| g helper | Primary, secondary, and tertiary amine definition | [1g/amine_descriptor.py](1g/amine_descriptor.py) |

These are copies of the generating routines with the input and output paths adapted for direct execution from each panel folder. They include the calculations, plot construction, colors, axes, and validation. Each script uses bundled inputs and writes to `generated_from_data` beside its code. All locations resolve absolutely from the script file, so execution does not depend on the working directory. There are no command-line options.

The shared modules remain the source used by the existing whole-figure export pipeline. Changes to a panel-local copy do not silently alter that pipeline. The existing final panel SVG, PDF, PNG, and Illustrator files are preserved. Native typography and layout adjustments are reproduced by the existing artwork-export wrappers.

The remaining final panels were retained from the supplied native artwork in this revision. Their export wrappers remain in the corresponding directories. Historical scripts for earlier layouts are not relabeled as generating code for the current final panels.

All three added plotting scripts were executed successfully from `/private/tmp`. Checks cover the unique hit counts, every inhibition density and width, the four property groups, all 125 structures' descriptors, 528 plotted points, axis limits, and the 33 amine examples. `figure_1_panel_source_code_manifest.json` records source provenance and file checksums.
