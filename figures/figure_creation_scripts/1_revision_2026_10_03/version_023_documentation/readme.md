# Figure 1 revision methods, version 023

The current property comparison contains log10 IPC, logP, and aromatic ring count for all hits in WT, Δ40, Δ3, and Δ6. The groups contain 8, 5, 39, and 124 compounds and overlap. Inactives are not plotted. Every compound is shown consistently as a point in every group and property, including observations beyond the box-plot whiskers. No plotted values are clipped.

## Current reproduction sequence

1. Run `isolated_panel_g_2026_10_04/methods/reproduce_isolated_properties.py` and `validate_isolated_properties.py` to regenerate and independently check the standalone SVGs and their source data. These paths are relative to the revision directory for documentation purposes. Executable scripts use absolute paths and require no command-line arguments.
2. Run `methods/prepare_three_property_panel.py` to lay out those plots in the existing panel g region.
3. Run `methods/integrate_three_properties_in_illustrator.jsx` in Illustrator to replace panel g in version 022 and save `Figure1_v023_three_properties.ai`. The script preserves other edits in an already-open source document and is intended to run once per source document.
4. Run `methods/export_three_property_figure.py` to generate the matching vector PDF, PNG and 300 dpi review PDF and verify that the artwork outside panel g matches the saved version 022.

Use `/opt/anaconda3/envs/Chemprop_Environment/bin/python` for analysis and plotting. Use `/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/venv/bin/python` for the PDF export script. All outputs go to `/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision`.

`preview_three_property_integration.py` creates a placement preview from the saved version 022 PDF artwork. This preview does not update the native Illustrator document. The final exports must come from the saved version 023 native file.

## Retained figure panels

The updated unique-hit counts, original replicate panels, and chemical-space maps are retained.

## Historical code

`reproduce_revision.py`, `violin_panels.py`, and `property_boxplots.py` reproduce earlier figure revisions with six property panels and different comparison groups. Run the current reproduction sequence above for the three-property comparison. Historical statistics are not analyses of the overlapping strain-hit groups in version 023. Earlier native-update scripts and descriptor definitions remain available for provenance. Documentation for version 022 is archived in `methods/version_022_documentation`.

## Illustrator save checks

The IPC subscript uses Arial in the generated panel PDF and in Illustrator. The integration workflow includes a font safeguard and suppresses interaction dialogs during the save. This avoids missing subscript glyphs and app-activation waits.

The native save in Illustrator 30.8.2 introduces a small colour conversion in panel f. Its vector geometry and opacity are identical to version 022, and the largest measured rendering difference is 4 levels on a 0–255 colour channel. Artwork above and below panel g is pixel-identical. The export script checks these bounds and records the full comparison in `tables/three_property_native_validation.json`.
