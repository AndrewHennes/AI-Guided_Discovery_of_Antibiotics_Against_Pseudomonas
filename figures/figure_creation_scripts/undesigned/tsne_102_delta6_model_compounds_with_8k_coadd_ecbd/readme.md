# Joint t-SNE of 102 compounds with 8K, CO-ADD, and ECBD

The joint analysis includes every valid standardized structure in the filtered 8K library, the latest filtered CO-ADD and ECBD datasets with the selected compounds removed, and the new 102-compound screening set. The complete CO-ADD and ECBD source records are used without activity filtering. All library structures share one gray background layer. The 102 new compounds retain their Δ6 activity highlighting, with 53 hits in magenta and 49 compounds below 80% inhibition as open black circles. Their one missing WT measurement remains missing.

The final embedding contains 206,734 distinct standardized structures. The gray background comprises 206,637 structures; five of the 102 highlighted structures also occur in the background libraries. Six CO-ADD source records have unparseable SMILES and two have missing SMILES; all are listed in `analysis/excluded_source_records.csv`. All ECBD, 8K, and new screening records are represented.

## Sources

- `Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv`, 7,997 records.
- `COADD_Antibiotic_Activity_Dataset_Filtered_Selected_Compounds_Removed.csv`, 185,495 records.
- `Formatted_ECBD_Antibiotic_Activity_Dataset_Filtered_Selected_Compounds_Removed.csv`, 100,761 records.
- `PA01_PA∆6_Growth_Inhibition_∆6_Model_Only_Processed.csv`, 102 records.

`inputs` contains losslessly compressed complete source CSVs, SHA-256 checksums, and the earlier removal records identifying the CO-ADD and ECBD dataset versions.

## Reproduce the figure

Run `plot_standardized_102_compounds_with_8k_coadd_ecbd_tsne.py` to recreate the SVG, PDF, and PNG using the saved joint coordinates. No command-line arguments are needed. Scripts resolve every input and output from their own absolute location and can run from any working directory.

Run `prepare_standardized_tsne.py` to repeat molecular standardization, fingerprint generation, principal component analysis, and t-SNE fitting. Alternatively, set `RECOMPUTE_EMBEDDING = True` in `config.py` before running the plotting script. `DATASET_PATHS` remains blank by default; the bundled source snapshots are used. External dataset paths, if supplied, must be absolute. Package versions are pinned in `requirements.txt`.

## Molecular preparation

All four sources use their original `SMILES` column and the same RDKit standardization procedure as the preceding 8K-only background analysis. This consists of sanitized parsing, `rdMolStandardize.Cleanup`, parent-fragment selection, `Uncharger`, explicit-hydrogen removal, and canonical isomeric SMILES serialization. Stereochemistry is retained. Cleanup and parent preparation repeat to a fixed point; charge-redistribution cycles use the lexicographically first canonical SMILES in the cycle. Tautomer canonicalization is not applied.

The standardization cache is keyed by the source checksums, procedure-code checksum, and RDKit version. `analysis/excluded_source_records.csv` identifies source records without a usable standardized structure. All other source records map to a joint embedding coordinate. Each standardized structure occurs once in the embedding, and dataset memberships preserve cross-library overlap. A structure shared with the new screening set receives its highlighted marker over the gray background.

## Embedding and graphical settings

Morgan binary fingerprints use radius 2 and 2,048 bits, with chirality excluded from the fingerprint as in the prior plots. Principal component analysis uses 50 components, the randomized solver, and seed 42. openTSNE uses perplexity 30, PCA initialization, Euclidean distance, Annoy neighbors, FFT optimization, 250 early-exaggeration iterations, 500 subsequent iterations, seed 42, and eight worker threads. All four datasets are fitted jointly.

Canvas dimensions, marker radii, outlines, gray and magenta colors, transparency, Helvetica typography, and the borderless layout follow Figure 1g and the preceding 102-compound plot. Coordinates are centered and scaled isotropically. No CO-ADD, ECBD, or 8K activity values affect background colors or membership. New-compound hit calls use inhibition ≥80%; missing measurements remain missing.

## Analysis products

- `102_compounds_standardized_with_tsne_coordinates.csv` retains the original 102 records with standardized SMILES, binary hit calls, and coordinates.
- `analysis/*_standardized_records.csv.gz` preserves the source row indices, identifiers, original SMILES, standardized structures, and processing status. The new screening set retains all original measurements.
- `analysis/joint_tsne_coordinates.csv` contains one row per standardized structure, dataset membership, and joint t-SNE coordinates.
- `analysis/*_records_with_tsne_coordinates.csv.gz` maps source records back to joint coordinates.
- `analysis/morgan_fingerprints.npz` stores fingerprint rows in `analysis/unique_standardized_molecules.csv` order.
- `analysis/pca_projection.npz` stores reduced coordinates, PCA components, mean, and explained variance ratios in that same order.
- `analysis/standardization_summary.json`, `analysis/analysis_settings.json`, `plot_validation.json`, and the export validation record counts, settings, and checks.
- `style_reference` contains the prior Figure 1g PDF and PNG used for graphical comparison.

SVG labels are vector outlines, the PDF embeds fonts, and plotted markers remain vector components.

The SVG export helper in `svg_export/` removes unused artwork and resolves all clipping and SVG masks into the exported paths and image data.
