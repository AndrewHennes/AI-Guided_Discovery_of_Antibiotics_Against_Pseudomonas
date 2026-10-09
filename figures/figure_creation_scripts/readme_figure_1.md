# Current Figure 1

This workflow reproduces the version 025 figure, including all active compounds in each strain group. It supersedes the earlier a–g layout. The current figure contains panels a–i.

## Reproduce

Run `reproduce_figure_1.py` from this folder using Python with the packages in `requirements.txt`. The script can be called by absolute path from any working directory and has no command-line arguments. An individual script in `1a` through `1i` refreshes its corresponding panel. Running the main script refreshes the complete figure and all publication SVGs.

Each panel folder receives SVG, PDF, and PNG files. `full_figure` receives the complete figure in those formats. Matching SVG copies are automatically written to `../figure_svg_files`. `current_figure_1_data` contains the regenerated data tables, direct plot exports, descriptor audits, and validation reports. The native Illustrator document is not edited by this workflow.

## Panel map

| Panel | Publication SVG | Entry script |
| --- | --- | --- |
| a | Figure 1a (Screening workflow).svg | `1a/antipseudomonal_compound_screening_workflow_a.py` |
| b | Figure 1b (Replicate correlations).svg | `1b/replicate_growth_inhibition_correlations_across_four_strains_b.py` |
| c | Figure 1c (Replicate correlation insets).svg | `1c/replicate_growth_inhibition_zoomed_views_c.py` |
| d | Figure 1d (Dataset hit counts).svg | `1d/unique_hit_counts_across_datasets_d.py` |
| e | Figure 1e (Hit overlap).svg | `1e/four_strain_antibacterial_hit_overlap_venn_e.py` |
| f | Figure 1f (Growth inhibition distributions).svg | `1f/growth_inhibition_log_density_violins_f.py` |
| g | Figure 1g (TPSA amines LogP).svg | `1g/TPSA_amine_count_LogP_all_strain_hits_g.py` |
| h | Figure 1h (8K COADD chemical space).svg | `1h/molecular_chemical_space_8k_coadd_tsne_h.py` |
| i | Figure 1i (8K Broad39K ECBD chemical space).svg | `1i/molecular_chemical_space_8k_broad39k_ecbd_tsne_i.py` |

## Current analysis and appearance

Panel d recalculates unique positive structures from the supplied filtered 8K, 39K, ECBD, and CO-ADD data. The 8K counts are WT 8, Δ40 5, Δ3 39, and Δ6 124. The other datasets contain 8 WT hits in 39K, 1 WT hit in ECBD, and 7 WT plus 75 Δ5 hits in CO-ADD. Duplicate records contribute once per canonical structure and strain.

Panel g contains **TPSA, number of amines, and LogP** as box-and-whisker plots. Each group includes all compounds active against that strain. Shared hits occur in every corresponding group. Counts are **WT 8, Δ40 5, Δ3 39, and Δ6 124**. All compounds are drawn with the same point style and deterministic horizontal jitter. There is no inactive group. Green, orange, red, and purple match the native figure. Y limits are −8 to 220 Å², −0.2 to 3.25 amines, and 0 to 8 LogP, respectively, and contain every observation.

Amine count combines primary, secondary, and tertiary amines, including protonated forms. The precise structural definition and 33 functional-group checks are in `shared/figure1_amine_descriptor.py`. Per-class counts are included in the output table.

Panel f retains the native version 025 artwork unchanged. The pipeline independently recomputes and checks all inhibition values and log-density curves against the verified reference tables.

Panels a, b, c, e, h, and i retain their native artwork. Panels h and i use the original coordinates, compound membership, and highlighted layers. The source data, embedding code, and saved coordinates remain in `shared/saved_plot_data` and the historical pipeline archive. Updated dataset counts refer to the current supplied files.

## Inputs and code

`figure1_current_config.py` defines current paths, panel boundaries, and export names. `shared/figure1_current_analysis.py` recalculates counts and inhibition densities; `shared/figure1_properties.py` calculates descriptors and draws panel g; `shared/figure1_current_rendering.py` assembles the figure and exports the standalone panels.

`shared/figure1_current/inputs` contains compressed tabular snapshots of the four supplied CSV files with a source manifest. No download or external chat output is required. `reference_tables` contains the verified counts, density curves, and descriptor table. `figure1_v025_reference.pdf` supplies the native vector layout. The analysis checks the supplied data against the verified reference values and stops if they disagree. Replacing the inputs for a new dataset requires reviewing and updating the reference expectations, rather than silently relabeling this figure as version 025.

The SVG exporter outlines text, embeds required imagery, removes off-canvas and fully covered graphics, and resolves SVG masks and clipping into geometry where supported. Complex native density curves retain standard SVG clipping when geometric intersection cannot be performed reliably. Each panel is standalone, with no external resource links. The standalone PDFs are generated from the cleaned vector SVGs.

## Validation

`figure_1_current_validation.json` records exact pixel comparisons for the retained panel f artwork, left panels, and chemical-space maps, plus the SVG structure checks and publication-file checksums. `current_figure_1_data/analysis_validation.json` records source checksums and the current counts. Panel g additionally checks membership independently against source-positive structure sets, all 125 unique structures' descriptor values, all 528 plotted points, axis limits, and the amine examples. Panel f checks every recalculated density and displayed width against its verified table.
