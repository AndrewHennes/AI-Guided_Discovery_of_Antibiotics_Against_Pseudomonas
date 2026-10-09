"""Plot standardized model-selected compounds with the filtered 8K library."""

import json
import os
from pathlib import Path
from svg_export import finalize_svg

import numpy as np
import pandas as pd

import config

DIRECTORY = config.DIRECTORY
os.environ.setdefault("MPLCONFIGDIR", str(DIRECTORY / "analysis" / "matplotlib_cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


def main():
    coordinate_path = DIRECTORY / "analysis" / "joint_tsne_coordinates.csv"
    if config.RECOMPUTE_EMBEDDING or not coordinate_path.is_file():
        from prepare_standardized_tsne import main as prepare
        prepare()
    frame = pd.read_csv(coordinate_path, float_precision="round_trip")
    style = config.PLOT_STYLE
    xy = frame[["tSNE_1", "tSNE_2"]].to_numpy()
    bounds = np.asarray(style["plot_bounds_pt"], dtype=float)
    center = (xy.min(axis=0) + xy.max(axis=0)) / 2
    scale = min(np.diff(bounds, axis=1).ravel() / np.ptp(xy, axis=0))
    display_xy = (xy - center) * scale + bounds.mean(axis=1)
    display_xy[:, 1] = bounds[1].sum() - display_xy[:, 1]
    frame[["Display_X_pt", "Display_Y_pt"]] = display_xy
    new = frame.in_New102.astype(bool).to_numpy()
    background = frame.in_8K.astype(bool).to_numpy()
    wt = new & frame.WT_New102.eq(1).to_numpy()
    delta6 = new & frame.Delta6_New102.eq(1).to_numpy() & ~wt
    other = new & ~(wt | delta6)
    width, height = style["canvas_width_pt"], style["canvas_height_pt"]
    plt.rcParams.update({
        "font.family": style["font_family"], "font.size": style["legend_font_size_pt"],
        "text.color": style["text_color"], "pdf.fonttype": 42,
        "svg.fonttype": "path", "svg.hashsalt": "standardized-102-compounds-with-8k",
        "savefig.facecolor": "white", "axes.unicode_minus": False,
    })
    figure = plt.figure(figsize=(width / 72, height / 72))
    axes = figure.add_axes([0, 0, 1, 1])
    axes.set(xlim=(0, width), ylim=(height, 0))
    axes.set_axis_off()
    layers = [
        ("8K filtered library", background, style["background_color"], "none", style["background_radius_pt"], style["background_opacity"], 0),
        ("New102 below Delta6 cutoff", other, "white", style["outline_color"], style["foreground_radius_pt"], 1, style["foreground_outline_width_pt"]),
        ("New102 Delta6 hit", delta6, style["delta6_color"], style["outline_color"], style["foreground_radius_pt"], style["foreground_opacity"], style["foreground_outline_width_pt"]),
        ("New102 WT hit", wt, style["wt_color"], style["outline_color"], style["foreground_radius_pt"], 1, style["foreground_outline_width_pt"]),
    ]
    plotted_layers = []
    for order, (label, mask, fill, edge, radius, opacity, linewidth) in enumerate(layers):
        if mask.any():
            artist = axes.scatter(display_xy[mask, 0], display_xy[mask, 1], s=(2 * radius) ** 2, marker="o", facecolors=fill, edgecolors=edge, linewidths=linewidth, alpha=opacity, zorder=order + 1, rasterized=False)
            artist.set_gid(label.replace(" ", "_"))
        plotted_layers.append({"layer": label, "structures": int(mask.sum()), "radius_pt": radius, "fill": fill, "edge": edge, "opacity": opacity, "outline_width_pt": linewidth})
    legend = [
        ("8K filtered library", style["background_color"], "none", 1.15, style["background_opacity"]),
        ("Δ6 active · 102-compound set", style["delta6_color"], style["outline_color"], 2.16, style["foreground_opacity"]),
        ("Δ6 <80% · 102-compound set", "white", style["outline_color"], 2.16, 1),
    ]
    if wt.any():
        legend.insert(1, ("WT active · 102-compound set", style["wt_color"], style["outline_color"], 2.16, 1))
    first_y = 294 - (len(legend) - 1) * 14.5
    for index, (label, fill, edge, radius, alpha) in enumerate(legend):
        y = first_y + index * 14.5
        axes.scatter([38], [y - 2.4], s=(2 * radius) ** 2, facecolors=fill, edgecolors=edge, linewidths=0 if edge == "none" else style["foreground_outline_width_pt"], alpha=alpha, zorder=10)
        axes.text(43, y, label, fontsize=style["legend_font_size_pt"], ha="left", va="baseline")
    for extension in ("svg", "pdf", "png"):
        options = {"dpi": 300} if extension == "png" else {}
        if extension == "svg":
            options["metadata"] = {"Date": None, "Description": "Joint t-SNE of standardized filtered 8K structures and 102 model-selected compounds"}
        figure.savefig(DIRECTORY / f"{config.FIGURE_STEM}.{extension}", **options)
        if extension == "svg":
            finalize_svg(DIRECTORY / f"{config.FIGURE_STEM}.svg")
    plt.close(figure)
    frame.to_csv(DIRECTORY / "analysis" / "plotted_coordinates.csv", index=False)
    report = {
        "reference": "Figure 1g molecular chemical-space panel",
        "style": style, "font_file_used": font_manager.findfont(style["font_family"], fallback_to_default=False),
        "layers": plotted_layers, "new_unique_structures_plotted": int(new.sum()),
        "new_record_count": sum(len(json.loads(value)) for value in frame.loc[new, "New102_Source_Rows_JSON"]),
        "coordinate_rescaling": {"center": center.tolist(), "scale": float(scale), "isotropic": True},
        "wt_hits_in_new_set": int(wt.sum()), "delta6_hits_in_new_set": int(delta6.sum()),
        "missing_wt_calls_in_new_set": int(frame.loc[new, "WT_New102"].isna().sum()),
        "legend": [row[0] for row in legend], "all_markers_vector": True,
    }
    (DIRECTORY / "plot_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {config.FIGURE_STEM}: SVG, PDF, PNG", flush=True)


if __name__ == "__main__":
    main()
