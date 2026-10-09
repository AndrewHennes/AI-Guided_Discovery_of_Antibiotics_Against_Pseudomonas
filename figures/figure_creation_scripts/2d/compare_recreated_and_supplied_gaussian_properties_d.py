"""Compare the 8K Gaussian fits with the supplied efflux-rule parameters."""

from io import BytesIO
from pathlib import Path
import json
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pymupdf
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject
from scipy.stats import norm

METHODS_DIRECTORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(METHODS_DIRECTORY))

import figure2_config as config
from shared.figure2_rendering import recolor, update_properties
from shared.rendering import compact_svg

OUTPUT_DIRECTORY = Path(__file__).resolve().parent / "gaussian_comparison"
PARAMETER_CATEGORIES = ("WT_and_Delta6_active", "Delta6_only_active")
PROPERTIES = ("log10IPC", "logP", "Number_of_Aromatic_Rings")
MEANS = np.array([[4.14, 5.89], [2.43, 3.75], [1.63, 2.55]])
STANDARD_DEVIATIONS = np.array([[2.07, 1.36], [1.22, 1.4], [1.4, 0.98]])
PRIORS = (0.25, 0.75)
PERCENTILES = (10, 25, 50, 75, 90)
PANEL_RECTANGLE = pymupdf.Rect(651, 19, 841.89, 207)


def supplied_fits():
    rows = []
    for i, prop in enumerate(PROPERTIES):
        for j, category in enumerate(PARAMETER_CATEGORIES):
            mean, sd = MEANS[i, j], STANDARD_DEVIATIONS[i, j]
            values = norm.ppf(np.asarray(PERCENTILES) / 100, loc=mean, scale=sd)
            rows.append({
                "Property": prop,
                "Hit_Category": category,
                "Gaussian_Mean": mean,
                "Gaussian_Sigma_MLE": sd,
                "Prior": PRIORS[j],
                **dict(zip([f"P{p}" for p in PERCENTILES], values)),
            })
    return pd.DataFrame(rows)


def comparison_table(current, supplied):
    rows = []
    for prop in PROPERTIES:
        for category in PARAMETER_CATEGORIES:
            a = current.loc[current.Property.eq(prop) & current.Hit_Category.eq(category)].iloc[0]
            b = supplied.loc[supplied.Property.eq(prop) & supplied.Hit_Category.eq(category)].iloc[0]
            row = {
                "Property": prop,
                "Hit_Category": category,
                "Current_N": int(a.N),
                "Current_Mean": a.Gaussian_Mean,
                "Supplied_Mean": b.Gaussian_Mean,
                "Mean_Change": b.Gaussian_Mean - a.Gaussian_Mean,
                "Current_SD": a.Gaussian_Sigma_MLE,
                "Supplied_SD": b.Gaussian_Sigma_MLE,
                "SD_Change_Percent": 100 * (b.Gaussian_Sigma_MLE / a.Gaussian_Sigma_MLE - 1),
            }
            for p in PERCENTILES:
                row[f"Current_P{p}"] = a[f"P{p}"]
                row[f"Supplied_P{p}"] = b[f"P{p}"]
                row[f"P{p}_Change"] = b[f"P{p}"] - a[f"P{p}"]
            rows.append(row)
    return pd.DataFrame(rows)


def adjust_common_axes(content):
    shift = (567.4322 - 440.1779) / 8 * 0.5
    baseline = 440.1779 - shift
    horizontal = b"q 1 0 0 1 682.1096 440.1779 cm\n0 0 m\n127.254 0 l"
    vertical = b"q 1 0 0 1 682.1096 440.1779 cm\n0 0 m\n0 127.254 l"
    if content.count(horizontal) != 1 or content.count(vertical) != 1:
        raise ValueError("Expected one horizontal and one vertical property axis.")
    content = content.replace(horizontal, f"q 1 0 0 1 682.1096 {baseline:.8f} cm\n0 0 m\n127.254 0 l".encode())
    content = content.replace(vertical, f"q 1 0 0 1 682.1096 {baseline:.8f} cm\n0 0 m\n0 {127.254 + shift:.8f} l".encode())
    for x in ("701.4099", "745.6307", "789.8516"):
        content = content.replace(f"q 1 0 0 1 {x} 440.1779 cm\n0 0 m\n0 -4.772 l".encode(), f"q 1 0 0 1 {x} {baseline:.8f} cm\n0 0 m\n0 -4.772 l".encode())
    start = content.index(b"9.7208 0 0 9.7208 683.0376 424.6582 Tm")
    end = content.index(b"ET", start)
    block = content[start:end]
    block = re.sub(rb"([\d.]+) ([\d.]+) Tm", lambda match: match[1] + b" " + f"{float(match[2]) - shift:.8f}".encode() + b" Tm", block)
    content = content[:start] + block + content[end:]
    content = content.replace(b"662.4067 470.2559 Tm", f"662.4067 {470.2559 - shift / 2:.8f} Tm".encode())
    return content


def make_panel(fits):
    page = PdfReader(config.ASSET_PATH).pages[0]
    content = update_properties(page.get_contents().get_data(), fits)
    content = recolor(adjust_common_axes(content), is_main=True)
    for ref in page["/Resources"]["/XObject"].values():
        form = ref.get_object()
        if form.get("/Subtype") == "/Form":
            form.set_data(recolor(form.get_data()))
    stream = DecodedStreamObject()
    stream.set_data(b"q\n" + content + b"\nQ\n")
    page[NameObject("/Contents")] = stream.flate_encode()
    writer = PdfWriter()
    writer.add_page(page)
    buffer = BytesIO()
    writer.write(buffer)
    source = pymupdf.open(stream=buffer.getvalue(), filetype="pdf")
    source[0].draw_rect(pymupdf.Rect(650, 12, 662, 43), color=None, fill=(1, 1, 1), overlay=True)
    source[0].draw_rect(pymupdf.Rect(681, 0, 821, 25), color=None, fill=(1, 1, 1), overlay=True)
    source[0].draw_rect(pymupdf.Rect(650, 197, 775, 208), color=None, fill=(1, 1, 1), overlay=True)
    result = pymupdf.open()
    target = result.new_page(width=PANEL_RECTANGLE.width, height=PANEL_RECTANGLE.height)
    target.show_pdf_page(target.rect, source, 0, clip=PANEL_RECTANGLE)
    target.draw_rect(pymupdf.Rect(0, 178, 124, 189), color=None, fill=(1, 1, 1), overlay=True)
    source.close()
    return result


def export(document, name, white_regions, dpi=216):
    stem = OUTPUT_DIRECTORY / name
    document.save(stem.with_suffix(".pdf"), garbage=3, deflate=True)
    svg = ET.fromstring(compact_svg(document[0].get_svg_image(text_as_path=True)))
    for x, y, width, height in white_regions:
        ET.SubElement(svg, "{http://www.w3.org/2000/svg}rect", {"x": str(x), "y": str(y), "width": str(width), "height": str(height), "fill": "#ffffff"})
    stem.with_suffix(".svg").write_text(ET.tostring(svg, encoding="unicode", xml_declaration=True), encoding="utf-8")
    document[0].get_pixmap(dpi=dpi, alpha=False).save(stem.with_suffix(".png"))


def main():
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    current = pd.read_csv(Path(__file__).resolve().parent / "Gaussian_fit_parameters_and_percentiles.csv")
    supplied = supplied_fits()
    supplied.to_csv(OUTPUT_DIRECTORY / "supplied_gaussian_parameters_and_percentiles.csv", index=False)
    comparison = comparison_table(current, supplied)
    comparison.to_csv(OUTPUT_DIRECTORY / "current_vs_supplied_gaussian_comparison.csv", index=False)
    panels = [make_panel(current), make_panel(supplied)]
    for panel, name in zip(panels, ("current_8k_gaussian_percentile_plot_d", "supplied_efflux_rule_gaussian_percentile_plot_d")):
        export(panel, name, [(0, 178, 124, 12)])
    margin, gap, title_height = 12, 20, 23
    width, height = PANEL_RECTANGLE.width, PANEL_RECTANGLE.height
    combined = pymupdf.open()
    page = combined.new_page(width=2 * margin + 2 * width + gap, height=2 * margin + title_height + height)
    for i, (panel, title) in enumerate(zip(panels, ("Current 8K fits", "Supplied Gaussians"))):
        x = margin + i * (width + gap)
        page.insert_textbox(pymupdf.Rect(x, margin, x + width, margin + title_height), title, fontsize=11, fontname="hebo", align=1)
        page.show_pdf_page(pymupdf.Rect(x, margin + title_height, x + width, margin + title_height + height), panel, 0)
    export(combined, "current_vs_supplied_gaussian_percentile_comparison_d", [(margin + i * (width + gap), margin + title_height + 178, 124, 20) for i in range(2)])
    for document in [*panels, combined]:
        document.close()
    metadata = {
        "parameter_column_mapping": {"Hits": "WT_and_Delta6_active", "Nonhits": "Delta6_only_active"},
        "IPC_scale": "log10IPC",
        "box_percentiles": [25, 75],
        "median_percentile": 50,
        "whisker_percentiles": [10, 90],
        "axis_limits": [-0.5, 8],
        "supplied_class_priors": {"Hits": 0.25, "Nonhits": 0.75},
        "priors_used_for_within_class_percentiles": False,
    }
    (OUTPUT_DIRECTORY / "comparison_settings.json").write_text(json.dumps(metadata, indent=2))
    print(comparison[["Property", "Hit_Category", "Mean_Change", "SD_Change_Percent"]].to_string(index=False))
    print(f"Saved comparison plots and tables to {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
