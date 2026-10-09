"""Fit independent four-parameter concentration-response curves for plate rows."""

import hashlib
import json
import os
from importlib.metadata import version
from pathlib import Path

import config
from svg_export import finalize_svg

ROOT = config.ANALYSIS_DIRECTORY
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / ".cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import expit
from scipy.stats import t

READOUT_TITLES = {
    "OD600": r"OD$_{600}$",
    "Resazurin_30min": "Resazurin, 30 min",
    "Resazurin_60min": "Resazurin, 60 min",
}
READOUT_SLUGS = {
    "OD600": "od600",
    "Resazurin_30min": "resazurin_30min",
    "Resazurin_60min": "resazurin_60min",
}
STRAIN_LABELS = {"WT": "PAO1 (WT)", "Delta6": "PAO1 (Δ6)"}
DOSES = config.HIGHEST_CONCENTRATION_UM / config.DILUTION_FACTOR ** np.arange(12)
CURVE_DOSES = np.geomspace(DOSES.min(), DOSES.max(), 500)


def read_plates():
    records, manifest = [], {}
    for readout, filename in config.INPUT_FILENAMES.items():
        supplied = config.INPUT_PATHS[readout]
        path = Path(supplied).expanduser() if supplied else config.INPUT_DIRECTORY / filename
        if not path.is_absolute():
            raise ValueError("Input paths must be absolute.")
        contents = path.read_bytes()
        lines = contents.decode("utf-16").splitlines()
        header = next(i for i, line in enumerate(lines) if "Temperature" in line)
        plate = []
        for line in lines[header + 1:]:
            cells = line.split("\t")
            if len(cells) < 14:
                continue
            try:
                values = [float(cell) for cell in cells[2:14]]
            except ValueError:
                continue
            plate.append(values)
        if np.asarray(plate).shape != (8, 12):
            raise ValueError(f"Expected eight rows and twelve columns in {filename}.")
        for row_index, values in enumerate(plate):
            row = "ABCDEFGH"[row_index]
            for col_index, raw in enumerate(values):
                records.append({
                    "Readout": readout,
                    "Strain": "WT" if row_index < 4 else "Delta6",
                    "Replicate": row_index % 4 + 1,
                    "Plate_Row": row,
                    "Plate_Column": col_index + 1,
                    "Well": f"{row}{col_index + 1:02d}",
                    "Concentration_uM": float(DOSES[col_index]),
                    "Raw_Signal": raw,
                    "Source_Filename": filename,
                })
        manifest[readout] = {
            "filename": filename,
            "sha256": hashlib.sha256(contents).hexdigest(),
            "plate_rows": 8, "plate_columns": 12,
            "encoding": "UTF-16LE with BOM",
        }
    (config.INPUT_DIRECTORY / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return pd.DataFrame(records)


def standardized_response(parameters, log_concentration):
    lower, log_amplitude, log_ic50, log_hill = parameters
    fraction = expit(np.exp(log_hill) * np.log(10) * (log_ic50 - log_concentration))
    return lower + np.exp(log_amplitude) * fraction


def fit_row(group):
    concentration = group.Concentration_uM.to_numpy()
    response = group.Raw_Signal.to_numpy()
    offset, scale = response.min(), np.ptp(response)
    if not np.isfinite(response).all() or scale <= 0:
        raise ValueError("Every row must have twelve finite, varying measurements.")
    standardized = (response - offset) / scale
    log_concentration = np.log10(concentration)

    def residuals(parameters):
        return standardized_response(parameters, log_concentration) - standardized

    lower_bounds = [-100, -10, log_concentration.min() - 6, np.log(0.01)]
    upper_bounds = [100, 10, log_concentration.max() + 6, np.log(100)]
    candidates = []
    for initial_ic50 in (5, 15, 50, 150):
        for initial_hill in (0.7, 1.5, 3):
            initial = [-0.05, np.log(1.1), np.log10(initial_ic50), np.log(initial_hill)]
            fit = least_squares(
                residuals, initial, bounds=(lower_bounds, upper_bounds),
                max_nfev=4000, ftol=1e-12, xtol=1e-12, gtol=1e-12,
            )
            if fit.success:
                candidates.append(fit)
    if not candidates:
        raise RuntimeError("No starting point converged.")
    fit = min(candidates, key=lambda candidate: candidate.cost)
    lower, log_amplitude, log_ic50, log_hill = fit.x
    df = len(response) - len(fit.x)
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * (2 * fit.cost / df)
    log_ic50_se = np.sqrt(max(0, covariance[2, 2]))
    ci_half_width = t.ppf(0.975, df) * log_ic50_se
    bottom = float(offset + scale * lower)
    top = float(offset + scale * (lower + np.exp(log_amplitude)))
    fitted = offset + scale * standardized_response(fit.x, log_concentration)
    near_best = [candidate.x[2] for candidate in candidates
                 if candidate.cost <= fit.cost + max(1e-9, fit.cost * 1e-6)]
    record = {
        "Readout": group.Readout.iloc[0], "Strain": group.Strain.iloc[0],
        "Replicate": int(group.Replicate.iloc[0]), "Plate_Row": group.Plate_Row.iloc[0],
        "IC50_uM": float(10 ** log_ic50),
        "IC50_95CI_Lower_uM": float(10 ** (log_ic50 - ci_half_width)),
        "IC50_95CI_Upper_uM": float(10 ** (log_ic50 + ci_half_width)),
        "Lower_Plateau_Raw": bottom, "Upper_Plateau_Raw": top,
        "Hill_Slope": float(np.exp(log_hill)),
        "R_squared": float(1 - np.sum((response - fitted) ** 2) / np.sum((response - response.mean()) ** 2)),
        "RMSE_Raw": float(np.sqrt(np.mean((response - fitted) ** 2))),
        "N_Concentrations": len(response), "Residual_Degrees_Of_Freedom": df,
        "Converged": bool(fit.success),
        "Optimizer_Bound_Active": bool(np.any(fit.active_mask)),
        "Near_Best_Starts_Log10_IC50_Range": float(np.ptp(near_best)),
    }
    measured = group.copy()
    measured["Fitted_Raw_Signal"] = fitted
    measured["Residual_Raw"] = response - fitted
    measured["Normalized_Inhibition_Percent"] = 100 * (top - response) / (top - bottom)
    measured["Fitted_Inhibition_Percent"] = 100 * (top - fitted) / (top - bottom)
    curve_raw = offset + scale * standardized_response(fit.x, np.log10(CURVE_DOSES))
    curve = pd.DataFrame({
        "Readout": record["Readout"], "Strain": record["Strain"],
        "Replicate": record["Replicate"], "Plate_Row": record["Plate_Row"],
        "Concentration_uM": CURVE_DOSES, "Fitted_Raw_Signal": curve_raw,
        "Fitted_Inhibition_Percent": 100 * (top - curve_raw) / (top - bottom),
    })
    return record, measured, curve


def configure_style():
    if Path(config.FONT_PATH).is_file():
        font_manager.fontManager.addfont(config.FONT_PATH)
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "axes.linewidth": 0.75, "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.major.size": 2.8, "ytick.major.size": 2.8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "svg.fonttype": "path", "svg.hashsalt": "princeton-bio-row-ic50",
        "pdf.fonttype": 42, "mathtext.fontset": "dejavusans",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def format_axis(ax, all_doses=True):
    ax.set_xscale("log")
    ax.set_xlim(DOSES.min() / 1.16, DOSES.max() * 1.16)
    ticks = DOSES[::-1] if all_doses else DOSES[::-1][::2]
    ax.set_xticks(ticks, [f"{dose:.3g}" for dose in ticks], rotation=35, ha="right")
    ax.minorticks_off()
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Concentration (µM)")


def draw_summary(ax, readout, measurements, curves, summary):
    subset = measurements[measurements.Readout == readout]
    curve_subset = curves[curves.Readout == readout]
    ax.axhline(50, color="#bcbcbc", lw=0.85, ls=(0, (5, 3)), zorder=0)
    handles, labels, low_values, high_values = [], [], [], []
    for strain in ("WT", "Delta6"):
        color, marker = config.COLORS[strain], config.MARKERS[strain]
        observed = subset[subset.Strain == strain].groupby("Concentration_uM").Normalized_Inhibition_Percent.agg(["mean", "std"])
        fitted = curve_subset[curve_subset.Strain == strain].groupby("Concentration_uM").Fitted_Inhibition_Percent.mean()
        stats = summary[(summary.Readout == readout) & (summary.Strain == strain)].iloc[0]
        ax.plot(fitted.index, fitted, color=color, lw=1.0)
        ax.errorbar(observed.index, observed["mean"], yerr=observed["std"],
                    fmt=marker, markersize=3.2, markerfacecolor=color,
                    markeredgecolor="#333333", markeredgewidth=0.5,
                    ecolor=color, elinewidth=0.8, capsize=1.5, capthick=0.65, zorder=3)
        ax.axvline(stats.IC50_Mean_uM, color=color, lw=0.75, ls=(0, (3, 5)))
        handles.append(Line2D([], [], marker=marker, linestyle="none", markersize=3.7,
                              markerfacecolor=color, markeredgecolor="#333333", markeredgewidth=0.5))
        labels.append(STRAIN_LABELS[strain] + "\n" + rf"IC$_{{50}}$ = {stats.IC50_Mean_uM:.2f} ± {stats.IC50_SD_uM:.2f} µM")
        low_values.extend((observed["mean"] - observed["std"]).tolist())
        high_values.extend((observed["mean"] + observed["std"]).tolist())
    ax.set_ylim(min(-10, min(low_values) - 5), max(110, max(high_values) + 5))
    ax.set_yticks(np.arange(0, 101, 20))
    format_axis(ax)
    ax.set_ylabel("Growth inhibition (%)")
    ax.set_title(READOUT_TITLES[readout], pad=8)
    legend = ax.legend(handles, labels, loc="upper left", frameon=False, fontsize=8,
                       handlelength=0.9, handletextpad=0.35, borderpad=0.1, labelspacing=0.55)
    for text, strain in zip(legend.get_texts(), ("WT", "Delta6")):
        text.set_color(config.COLORS[strain])


def save_figure(fig, stem):
    for extension in ("svg", "pdf", "png"):
        kwargs = {"dpi": config.PNG_DPI} if extension == "png" else {}
        if extension == "svg":
            kwargs["metadata"] = {"Date": None}
        fig.savefig(ROOT / f"{stem}.{extension}", **kwargs)
        if extension == "svg":
            finalize_svg(ROOT / f"{stem}.svg")
    plt.close(fig)


def make_plots(measurements, curves, row_fits, summary):
    configure_style()
    readouts = list(config.INPUT_FILENAMES)
    for readout in readouts:
        fig, ax = plt.subplots(figsize=(3.65, 3.25))
        fig.subplots_adjust(left=0.19, right=0.98, bottom=0.24, top=0.88)
        draw_summary(ax, readout, measurements, curves, summary)
        save_figure(fig, f"princeton_bio_{READOUT_SLUGS[readout]}_ic50")

        fig, axes = plt.subplots(4, 2, figsize=(7.4, 9.6))
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.065, top=0.93, hspace=0.69, wspace=0.35)
        fig.suptitle(READOUT_TITLES[readout], fontsize=12, y=0.98)
        for index, row in enumerate("ABCDEFGH"):
            row_data = measurements[(measurements.Readout == readout) & (measurements.Plate_Row == row)]
            row_curve = curves[(curves.Readout == readout) & (curves.Plate_Row == row)]
            fit = row_fits[(row_fits.Readout == readout) & (row_fits.Plate_Row == row)].iloc[0]
            ax = axes[index % 4, index // 4]
            color = config.COLORS[fit.Strain]
            ax.plot(row_curve.Concentration_uM, row_curve.Fitted_Raw_Signal, color=color, lw=1.0)
            ax.scatter(row_data.Concentration_uM, row_data.Raw_Signal, s=15,
                       marker=config.MARKERS[fit.Strain], facecolor=color, edgecolor="#333333", lw=0.5, zorder=3)
            ax.axhline(fit.Lower_Plateau_Raw, color="#999999", ls=(0, (2, 3)), lw=0.6)
            ax.axhline(fit.Upper_Plateau_Raw, color="#999999", ls=(0, (2, 3)), lw=0.6)
            ax.axvline(fit.IC50_uM, color=color, ls=(0, (3, 5)), lw=0.7)
            format_axis(ax, all_doses=False)
            ax.set_ylabel(r"OD$_{600}$" if readout == "OD600" else "Fluorescence (RFU)", fontsize=9)
            ax.set_title(f"{STRAIN_LABELS[fit.Strain]}, row {row}   " + rf"IC$_{{50}}$ = {fit.IC50_uM:.2f} µM", fontsize=8.5, color=color)
        save_figure(fig, f"princeton_bio_{READOUT_SLUGS[readout]}_individual_row_fits")

    fig, axes = plt.subplots(1, 3, figsize=(10.7, 3.25))
    fig.subplots_adjust(left=0.063, right=0.99, bottom=0.24, top=0.88, wspace=0.30)
    for ax, readout in zip(axes, readouts):
        draw_summary(ax, readout, measurements, curves, summary)
    save_figure(fig, "princeton_bio_three_readout_ic50_comparison")


def main():
    raw = read_plates()
    raw.to_csv(ROOT / "plate_measurements_long.csv", index=False)
    results = [fit_row(group) for _, group in raw.groupby(["Readout", "Strain", "Plate_Row"], sort=False)]
    row_fits = pd.DataFrame([item[0] for item in results])
    measurements = pd.concat([item[1] for item in results], ignore_index=True)
    curves = pd.concat([item[2] for item in results], ignore_index=True)
    summary = row_fits.groupby(["Readout", "Strain"], sort=False).IC50_uM.agg(
        N_Replicates="count", IC50_Mean_uM="mean", IC50_SD_uM="std", IC50_Median_uM="median"
    ).reset_index()
    row_fits.to_csv(ROOT / "ic50_per_row.csv", index=False)
    summary.to_csv(ROOT / "ic50_summary_by_readout_and_strain.csv", index=False)
    measurements.to_csv(ROOT / "row_normalized_inhibition_and_residuals.csv", index=False)
    curves.to_csv(ROOT / "fitted_curves_per_row.csv", index=False)
    make_plots(measurements, curves, row_fits, summary)
    validation = {
        "total_measurements": len(raw), "total_independent_row_fits": len(row_fits),
        "concentrations_uM_left_to_right": DOSES.tolist(),
        "all_fits_converged": bool(row_fits.Converged.all()),
        "any_optimizer_bound_active": bool(row_fits.Optimizer_Bound_Active.any()),
        "all_ic50_inside_tested_range": bool(row_fits.IC50_uM.between(DOSES.min(), DOSES.max()).all()),
        "minimum_row_r_squared": float(row_fits.R_squared.min()),
        "largest_near_best_start_log10_ic50_range": float(row_fits.Near_Best_Starts_Log10_IC50_Range.max()),
        "versions": {name: version(name) for name in ("numpy", "pandas", "scipy", "matplotlib")},
    }
    (ROOT / "fit_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
