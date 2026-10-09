# Figure 1 revision methods, version 024

The current panel g shows TPSA, number of amines, and LogP for 8 WT, 5 Δ40, 31 Δ3, and 85 Δ6 compounds. WT and Δ40 contain all active compounds in their respective strains. Δ3 excludes WT-active compounds. Δ6 excludes compounds active in either Δ3 or WT. These definitions appear in the caption, while plot labels use only the strain names and sample sizes.

## Reproduction

1. Run `reproduce_activity_group_panel.py` with `/opt/anaconda3/envs/Chemprop_Environment/bin/python` to regenerate the panel and validated data tables. `amine_descriptor.py` supplies the existing primary, secondary, and tertiary amine definition.
2. Run `integrate_activity_groups_in_illustrator.jsx` in Illustrator to replace only panel g in version 023, save version 024, and close the saved document. The script suppresses interaction dialogs during the edit and save.
3. Run `export_activity_group_figure.py` with `/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/venv/bin/python` to export the vector PDF and review PDF and verify that all artwork outside panel g, including panel f, is pixel-identical to version 023.

Scripts require no command-line arguments. All inputs and outputs use absolute paths. Results go to `/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision`. Data and panel files are in `panel_g_activity_groups_v024`. The original source file and previous native revisions are retained.

## Validation

The group definitions are checked independently against sets of positively labelled source records. Descriptor values are compared with the previously verified table, and 33 amine examples are checked. The plot audit verifies all 387 compound points and their axis coverage. Native export checks include the corrected labels, 12 count annotations, Arial fonts, and exact comparisons outside panel g. Panel f is not regenerated or edited.

The previous scripts reproduce earlier figure versions and cohort definitions. Documentation for version 023 is archived in `methods/version_023_documentation`. See `METHODS_AND_CAPTION.md` and the version 024 panel folder for the current definitions and caption.
