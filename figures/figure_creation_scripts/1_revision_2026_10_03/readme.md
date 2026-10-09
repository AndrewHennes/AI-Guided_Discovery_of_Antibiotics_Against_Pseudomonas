# Figure 1 revision methods, version 025

The current panel g shows TPSA, number of amines, and LogP for all active compounds in each strain. Counts are 8 WT, 5 Δ40, 39 Δ3, and 124 Δ6. Shared hits appear in every strain group where they are active. Plot labels use only strain names and counts; the caption explains the overlapping groups.

## Reproduction

1. Run `reproduce_all_hit_property_panel.py` with `/opt/anaconda3/envs/Chemprop_Environment/bin/python` to regenerate the panel and validated data tables. `amine_descriptor.py` supplies the existing primary, secondary, and tertiary amine definition.
2. Run `integrate_all_hits_in_illustrator.jsx` in Illustrator to replace only panel g and save version 025. It prefers the already-open Downloads copy of version 024 and otherwise uses the delivered version 024. Existing unsaved edits are carried forward into the new file. The script captures native PNG exports before and after editing for comparison.
3. Run `export_all_hit_figure.py` with `/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/venv/bin/python` to export the saved native figure and verify the labels, counts, fonts, and preservation of panel f and the other artwork.

Scripts use absolute input and output paths and require no command-line arguments. Outputs go to `/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision`. Data and panel artifacts are in `panel_g_all_hits_v025`.

## Validation

Membership is checked against sets of positively labelled source records, without cross-strain exclusions. The groups have 176 observations across 125 unique structures. The audit verifies all 528 points, descriptor values, and axis coverage. Native before/after images must match exactly outside panel g, including panel f. Comparisons with the previously delivered version 024 are also recorded, while the open document provides the authoritative baseline for preserving any user edits.

Earlier scripts remain available for historical revisions. Version 024 documentation is archived in `methods/version_024_documentation`.
