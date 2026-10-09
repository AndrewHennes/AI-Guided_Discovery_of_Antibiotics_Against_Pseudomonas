# Princeton Bio concentration-response analysis

The three supplied plates measure optical density at 600 nm (OD600) and resazurin fluorescence at 30 and 60 minutes. Rows A–D are four WT replicates; rows E–H are four Δ6 replicates. Column 1 is 255 µM, and each successive column halves the concentration through column 12 at 0.12451171875 µM.

Each row is fitted independently to a four-parameter logistic model using all twelve raw measurements and equal residual weights:

\[
y(c)=B+\frac{T-B}{1+(c/\mathrm{IC}_{50})^h}.
\]

The lower plateau \(B\), upper plateau \(T\), half-maximal inhibitory concentration \(\mathrm{IC}_{50}\), and Hill slope \(h\) are fitted separately for every row. The upper plateau exceeds the lower plateau, and the slope and IC50 are positive. Neither plateau is fixed to a measured well or to zero. The optimizer uses twelve initial combinations of IC50 and slope; the converged solution with the smallest squared error is retained.

The displayed inhibition for each measured well is calculated using its own row's fitted plateaus:

\[
\mathrm{Inhibition}(c)=100\frac{T-y(c)}{T-B}.
\]

Values are retained without clipping. Comparison-plot points and error bars show the arithmetic mean and sample standard deviation of the four normalized rows. Solid lines show the pointwise arithmetic mean of the four independently fitted normalized curves. Colored vertical dashed lines and legend values show the arithmetic mean of the four row-level IC50 values; the accompanying uncertainty is their sample standard deviation. The gray horizontal dashed line marks 50% inhibition.

Individual-row plots show the raw measurements and fitted curve for each row. Gray dotted horizontal lines mark its fitted plateaus, and the colored vertical line marks its IC50. Fluorescence is expressed in relative fluorescence units (RFU).

## Results

| Readout | WT IC50, µM (mean ± SD) | Δ6 IC50, µM (mean ± SD) |
| --- | ---: | ---: |
| OD600 | 19.55 ± 2.52 | 17.92 ± 5.09 |
| Resazurin, 30 min | 8.53 ± 1.32 | 5.88 ± 1.12 |
| Resazurin, 60 min | 19.32 ± 1.74 | 31.43 ± 4.20 |

`ic50_per_row.csv` contains all 24 IC50 values, fitted lower and upper plateaus, slopes, fit statistics, and approximate 95% confidence intervals calculated from the local covariance of log10(IC50) with eight residual degrees of freedom. `ic50_summary_by_readout_and_strain.csv` contains the four-replicate summaries.

## Files and reproduction

- `inputs/` contains the original three plate files and their SHA-256 checksums.
- `config.py` specifies the dilution series, strain colors, and filenames. External input paths are blank; the bundled input files are used by default.
- `fit_and_plot_princeton_bio_ic50.py` parses the plates, fits all 24 rows, writes all tables, and produces every SVG, PDF, and PNG in this folder. It resolves file locations from its own absolute location and runs from any working directory without command-line arguments.
- `plate_measurements_long.csv` contains the 288 raw measurements and their well, strain, replicate, and concentration assignments.
- `row_normalized_inhibition_and_residuals.csv` adds fitted values, residuals, and inhibition normalized using each row's fitted plateaus.
- `fitted_curves_per_row.csv` contains the 24 fitted curves sampled at 500 concentrations.
- `princeton_bio_three_readout_ic50_comparison.*` contains the three-readout comparison. Separate comparison plots are supplied for each readout.
- `princeton_bio_*_individual_row_fits.*` contains the individual raw-signal row fits for each readout.
- `fit_validation.json` records convergence, the concentration series, and software versions.
- `requirements.txt` lists the pinned Python packages.

WT uses green circles and Δ6 uses magenta squares. SVG labels are vector outlines and the PDF files embed fonts. Each SVG contains only its displayed plot or comparison, with no unrelated figure artwork.

The SVG export helper in `svg_export/` removes unused artwork and resolves all clipping and SVG masks into the exported paths and image data.
