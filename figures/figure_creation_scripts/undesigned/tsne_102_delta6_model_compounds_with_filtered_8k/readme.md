# Standardized 102-compound and filtered 8K t-SNE

This analysis jointly embeds the 102 compounds in `PA01_PA∆6_Growth_Inhibition_∆6_Model_Only_Processed.csv` and the 7,997 records in `Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv`. The joint embedding contains 8,092 distinct standardized structures, comprising 7,990 from 8K and all 102 new compounds. There is no standardized structure overlap between the two sets.

The 8K library is displayed as small gray points. The new compounds are displayed as larger circles, with Δ6 hits in magenta and Δ6 inhibition below 80% as open black circles. There are 53 Δ6 hits and 49 compounds below the Δ6 cutoff. Figure dimensions, background and foreground marker sizes, line widths, alpha values, typeface, and the borderless layout follow Figure 1g.

## Reproduction

Run `plot_standardized_102_compounds_with_filtered_8k_tsne.py` to regenerate the SVG, PDF, and PNG from the included t-SNE coordinates. No command-line arguments or external dataset paths are needed. All locations resolve from the script's absolute location, so it can run from any working directory.

Run `prepare_standardized_tsne.py` to repeat SMILES standardization, fingerprint generation, principal component analysis, and t-SNE fitting from the saved source files. Alternatively, set `RECOMPUTE_EMBEDDING = True` in `config.py` before running the plotting script. The package versions used are pinned in `requirements.txt`. `DATASET_PATHS` is blank by default; external paths, when supplied, must be absolute. Plot settings are in `PLOT_STYLE`.

## Molecular preparation and embedding

Both sources use their original `SMILES` column as input to the same procedure. RDKit sanitizes the structures, applies `rdMolStandardize.Cleanup`, selects the parent fragment, applies `Uncharger`, removes explicit hydrogens, and writes canonical isomeric SMILES. Stereochemistry is retained. Cleanup and parent preparation repeat to a fixed point; a charge-redistribution cycle is assigned the lexicographically first canonical SMILES in that cycle. Tautomer canonicalization is not applied. Every resulting structure passes a second-standardization consistency check.

Each distinct standardized structure receives one embedding coordinate. Per-record outputs preserve all source identifiers and measurements. The 8K binary WT and Δ6 columns are retained; new WT and Δ6 binary columns use growth inhibition ≥80%, preserving missing values. Repeated measurements of one structure aggregate by the maximum available binary call for each strain.

Morgan binary fingerprints use radius 2 and 2,048 bits, with chirality excluded from the fingerprint as in the earlier plots. Principal component analysis reduces the fingerprints to 50 components using the randomized solver and seed 42. openTSNE uses perplexity 30, PCA initialization, Euclidean distance, Annoy neighbors, FFT optimization, 250 early-exaggeration iterations, 500 subsequent iterations, seed 42, and eight worker threads. These settings match the earlier t-SNE generation code. Both datasets are fitted together. Plot coordinates are centered and scaled isotropically to the Figure 1g plotting extent.

## Files

- `inputs` contains losslessly compressed copies of both complete original CSVs and their checksums.
- `102_compounds_standardized_with_tsne_coordinates.csv` contains all new compounds, original measurements, standardized SMILES, binary calls, and coordinates.
- `analysis/*_standardized_records.csv.gz` contains all original source columns plus standardized structures and hit calls.
- `analysis/joint_tsne_coordinates.csv` contains the joint structure-level embedding and dataset membership.
- `analysis/*_records_with_tsne_coordinates.csv.gz` contains each source record joined to its coordinate.
- `analysis/morgan_fingerprints.npz` and `analysis/pca_projection.npz` contain fingerprints and PCA inputs and parameters, keyed in the joint molecule-table order.
- `analysis/standardization_changes.csv` lists every record changed during standardization.
- `analysis/standardization_summary.json`, `analysis/analysis_settings.json`, and `plot_validation.json` record counts, settings, and checks.
- `style_reference` contains the prior Figure 1g PDF and PNG used to match the graphical style.

SVG text is stored as vector outlines. The PDF embeds fonts. Markers remain individual vector components.

The SVG export helper in `svg_export/` removes unused artwork and resolves all clipping and SVG masks into the exported paths and image data.
