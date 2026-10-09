# Figure 1 generation-code map

The code from this task is organized according to the panel positions in `Figure1_v025_all_hits_TPSA_amines_LogP (3).ai`.

| Final panel | Plot | Generation code |
| --- | --- | --- |
| 1a | Screening workflow | `1a/source_code` |
| 1b | Replicate correlations | `1b/source_code` |
| 1c | Replicate correlation zoomed views | `1c/source_code` |
| 1d | Unique hits across four datasets | `1d/source_code` |
| 1e | Four-strain hit overlap | `1e/source_code` |
| 1h | 8K and CO-ADD chemical space | `1h/source_code` |
| 1i | 8K, Broad 39K, and ECBD chemical space | `1i/source_code` |

The existing main scripts in `1a` through `1i` export the final layout. The added `source_code` directories retain the underlying generation implementations and provide separate reproduction entry scripts. The current panel 1f inhibition-density and panel 1g physicochemical-property workflows are already implemented by `shared/figure1_inhibition.py`, `shared/figure1_current_analysis.py`, `shared/figure1_properties.py`, and `shared/figure1_amine_descriptor.py`.

The chemical-space source maps use the saved coordinates already included in `shared/saved_plot_data`. The source reproduction writes outputs into `/Users/asselism/Documents/Codex/2026-09-21/please-create-an-html-file-that/outputs/final_figure1_source_reproduction`. `figure1_generation_code_manifest.json` records the source files and checksums.
