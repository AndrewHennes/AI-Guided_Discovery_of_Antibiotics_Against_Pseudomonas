# Unique Pseudomonas hits by dataset

## October 3 update

The count tables and SVG/PDF/PNG count plots were regenerated from the supplied selected-compound-removed ECBD, CO-ADD, and 39K files, plus the supplied filtered 8K file. WT counts are now 8, 8, 1, and 7 for 8K, 39K, ECBD, and CO-ADD, respectively. CO-ADD retains 75 unique Δ5 hit structures. Purchased Broad counts are unchanged. The script reads absolute snapshot paths and writes products to the revision output directory in the Codex workspace. The source manifest is in `../1_revision_2026_10_03/saved_inputs/source_manifest.json`.

Use the refreshed CSV tables for the October 3 count analysis. The explorer and the following workflow use the earlier dataset snapshots.

## Historical workflow and explorer documentation

Run `plot_unique_pseudomonas_hit_counts.py` to regenerate the count tables and the SVG, PDF, and PNG plot. The script resolves every input and output from its absolute location and has no command-line arguments. Dataset paths are blank by default; the included tables in `saved_inputs` contain all source records and the columns needed for this comparison. To use external source files, enter absolute paths in `DATASET_PATHS`.

The plot uses the filtered Broad 8K, Broad 39K, ECBD, and CO-ADD datasets used in Figure 1, plus the 183 purchased Broad compounds from the supplied screening workbook. The Broad library prediction workbook contains model scores and does not supply additional measured hits.

A unique molecule is a canonical isomeric molecular structure within a dataset. Canonical SMILES are computed from `Standardized_SMILES` when available, otherwise `SMILES`, retaining stereochemistry and disconnected components. Each structure contributes once to each measured strain. A structure with any binary hit call equal to 1 is a hit. Existing dataset binary columns supply the 80% growth inhibition calls; purchased Broad WT and Δ6 calls use measured inhibition strictly greater than 80%. Missing calls remain missing.

| Dataset | WT | Δ40 | Δ3 | Δ6 | Δ5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Broad 8K | 8 | 5 | 39 | 124 | |
| Broad 39K | 17 | | | | |
| ECBD | 3 | | | | |
| CO-ADD | 14 | | | | 75 |
| Broad purchases | 5 | | | 53 | |

Mutant totals include molecules that are also WT hits. Blank entries indicate that the dataset does not include that Pseudomonas mutant assay. CO-ADD's Δ5 hits represent 75 unique structures and 76 compound IDs. Both counts are retained in the detailed CSV; the plot uses structures consistently across datasets. There are 41 distinct WT hit structures across the four filtered datasets and 46 after adding the purchased Broad compounds. One WT hit structure occurs in both 8K and 39K.

`WT_and_mutant_hit_counts.csv` is the concise count table. `unique_hit_counts_by_dataset_and_strain.csv` includes structure counts, source-record counts, compound-ID counts where available, assay columns, and measured structure counts. `unique_hit_compounds.csv` identifies every plotted hit structure and its source IDs and row indices. `WT_hit_overlap_between_datasets.csv` records cross-dataset overlap. Original filenames, file checksums, and input row counts are in `saved_inputs/source_manifest.json`. The validation report compares the computed counts against an independent read of the source datasets.

The plot uses the Figure 1 strain colors, Arial text, black bar outlines, and integer labels above each bar. SVG text is represented as vector outlines, and the PDF embeds editable text fonts.

`Pseudomonas_hit_explorer.html` is a self-contained offline explorer with chemical structures, dataset and strain filters, text search, a table view, molecular details, and export of the filtered results as CSV. It contains 286 dataset entries representing 285 distinct hit structures. The initial WT view contains 47 entries representing 46 distinct structures.

Run `build_hit_explorer.py` to regenerate the HTML and its embedded data. External dataset paths are blank. The included `saved_inputs/*_explorer_records.csv.gz` files contain every source record belonging to a hit structure, including associated nonhit and missing calls, identifiers, and inhibition measurements. Their provenance is recorded in `saved_inputs/explorer_source_manifest.json`. The script generates molecular drawings and properties with RDKit and combines the data with `hit_explorer_template.html`. Inputs and outputs resolve from the script's absolute location, independently of the current working directory.

The SVG export helper in `svg_export/` removes unused artwork and resolves all clipping and SVG masks into the exported paths and image data.
