"""Vector geometry and exports for Figure 2."""

from io import BytesIO
import re

import pymupdf
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject

from config import COLORS_CMYK
import figure2_config as config
from shared import figure2_analysis as analysis
from shared.rendering import compact_svg, numbers

COLOR_PATTERN = re.compile(rb"([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+) +([sS][cC][nN])")
OLD_MAGENTA = (0.357, 0.816, 0, 0)
WT_COLORS = {(1, 0, 1, 0), (0.832, 0.086, 0.887, 0.004), (0.816, 0.109, 0.855, 0.012), (0.813, 0.102, 0.875, 0.008), (0.644, 1, 0.195, 0.058)}
DELTA6_COLORS = {(0.457, 0.039, 0.105, 0), (0.457, 0.039, 0.121, 0), (0.625, 0.199, 0.137, 0), (0.527, 0, 0.125, 0), (0.637, 0, 0.207, 0), (0.102, 0.973, 0.98, 0.016)}


def recolor(content, is_main=False):
    def replace(match):
        old = tuple(float(v) for v in match.group(1, 2, 3, 4))
        target = None
        if old in WT_COLORS:
            target = COLORS_CMYK["WT"]
        elif old in DELTA6_COLORS:
            target = COLORS_CMYK["Delta6"]
        elif old == OLD_MAGENTA:
            following = content[match.end():match.end() + 160]
            keep = is_main and any(anchor in following for anchor in (b"32.8146 452.451", b"63.7597 465.2715", b"287.7709 152.5137", b"531.5226 319.9182"))
            if not keep:
                target = COLORS_CMYK["Delta6"]
        return numbers(target).encode() + b" " + match[5] if target is not None else match[0]
    return COLOR_PATTERN.sub(replace, content)


def curve_commands(curve):
    commands = ["q 1 0 0 1 32.8146 452.451 cm"]
    for i, (fpr, tpr) in enumerate(zip(curve["fpr"], curve["tpr"])):
        commands.append(numbers((115.508 * fpr, 115.508 * tpr)) + (" m" if i == 0 else " l"))
    commands.extend(("S", "Q"))
    return ("\n".join(commands) + "\n").encode()


def update_roc(content, curves):
    for model, color in (("GI60", b"0.883 0.43 0.418 0.121  SCN"), ("GI80", b"0.357 0.816 0 0  SCN")):
        if model == "GI60":
            color_start = content.index(color)
        else:
            color_start = content.index(color, content.index(b"0.883 0.43 0.418 0.121  SCN"))
        start = content.index(b"q 1 0 0 1 32.8146 452.451 cm\n", color_start)
        end = content.index(b"S\nQ\n", start) + len(b"S\nQ\n")
        content = content[:start] + curve_commands(curves[model]) + content[end:]
    return content


def property_commands(fits):
    baseline = 440.1779
    scale = (567.4322 - baseline) / 8
    lefts = {"Delta6_only_active": (687.9421, 732.1630, 776.4899), "WT_and_Delta6_active": (704.4852, 748.8121, 793.0330)}
    width = 10.392
    commands = ["q", "/CS0 cs /CS0 CS", "[] 0 d"]
    for category, color in (("Delta6_only_active", COLORS_CMYK["Delta6"]), ("WT_and_Delta6_active", COLORS_CMYK["WT"])):
        for prop, left in zip(analysis.PROPERTIES, lefts[category]):
            row = fits.loc[fits.Property.eq(prop) & fits.Hit_Category.eq(category)].iloc[0]
            p10, p25, p50, p75, p90 = [baseline + scale * row[p] for p in ("P10", "P25", "P50", "P75", "P90")]
            center = left + width / 2
            stroke_width = 0.75 if category == "Delta6_only_active" and prop == "Number_of_Aromatic_Rings" else 1
            commands.extend((numbers(color) + " scn", "0 0 0 1 SCN", f"{stroke_width} w 0 J [] 0 d", numbers((left, p75, width, p25 - p75)) + " re B"))
            commands.extend(("0 0 0 1 SCN", "0.75 w 0 J [] 0 d", numbers((left, p50)) + " m", numbers((left + width, p50)) + " l S"))
            commands.extend(("0 0 0 1 SCN", "1 w 1 J [] 0 d"))
            for y in (p10, p90):
                commands.extend((numbers((center - 2.7572, y)) + " m", numbers((center + 2.8628, y)) + " l S"))
            for start, end in ((p25 - 0.318, p10), (p90, p75 + 0.318)):
                commands.extend((numbers((center, start)) + " m", numbers((center, end)) + " l S"))
    commands.extend(("Q", "/CS0 CS 0 0 0 1 SCN 1 w 1 J [] 0 d"))
    return ("\n".join(commands) + "\n").encode()


def update_properties(content, fits):
    start = content.index(b"0.357 0.816 0 0  scn\n687.942")
    end = content.index(b"BT\n0.746 0.676 0.668 0.898  scn\n-0.016 Tc", start)
    content = content[:start] + property_commands(fits) + content[end:]
    old = b"9.7208 0 0 9.7208 693.3081 424.6582 Tm\n(IPC)Tj\n4.009 0 Td\n(cLogP)Tj\n3.633 -0.023 Td\n(Number of)Tj\n0.417 -1.2 Td\n(aromatic)Tj\n0.695 -1.2 Td\n(rings )Tj"
    new = b"9.7208 0 0 9.7208 683.0376 424.6582 Tm\n(log)Tj\n6.80456 0 0 6.80456 696.0076 422.9282 Tm\n(10)Tj\n9.7208 0 0 9.7208 703.5776 424.6582 Tm\n(IPC)Tj\n9.7208 0 0 9.7208 732.2787872 424.6582 Tm\n(cLogP)Tj\n3.633 -0.023 Td\n(Number of)Tj\n0.417 -1.2 Td\n(aromatic)Tj\n0.695 -1.2 Td\n(rings )Tj"
    if content.count(old) != 1:
        raise ValueError("Expected one property-axis text anchor.")
    return content.replace(old, new)


def assembled_pdf():
    curves = analysis.roc_tables()
    _, fits, comparisons = analysis.property_tables()
    if not comparisons.p_value.lt(0.05).all():
        raise ValueError("Update the significance annotations for these property comparisons.")
    page = PdfReader(config.ASSET_PATH).pages[0]
    content = page.get_contents().get_data()
    content = update_roc(content, curves)
    content = update_properties(content, fits)
    content = recolor(content, is_main=True)
    from shared.figure2_scatter_svg import remove_previous_scatter, add_scatter_components
    content = remove_previous_scatter(content)
    visited = set()

    def update_forms(resources):
        for name, reference in resources.get("/XObject", {}).items():
            form = reference.get_object()
            if id(form) in visited:
                continue
            visited.add(id(form))
            if form.get("/Subtype") != "/Form":
                continue
            data = form.get_data()
            if name in ("/Fm16", "/Fm17"):
                model = "GI60" if name == "/Fm16" else "GI80"
                data = re.sub(rb"AUC = 0\.\d{3}", f"AUC = {curves[model]['auc']:.3f}".encode(), data)
            form.set_data(recolor(data))
            if "/Resources" in form:
                update_forms(form["/Resources"])
    update_forms(page["/Resources"])
    stream = DecodedStreamObject()
    stream.set_data(content)
    page[NameObject("/Contents")] = stream.flate_encode()
    writer = PdfWriter()
    writer.add_page(page)
    writer.add_metadata({"/Title": "DeepFlux predictions and hit properties", "/Subject": "Delta5 hits versus inactive compounds using original ensemble scores; Gaussian property percentiles"})
    output = BytesIO()
    writer.write(output)
    return add_scatter_components(output.getvalue())


def export_artifacts(data, panel=None):
    directory = config.METHODS_DIRECTORY / f"2{panel}" if panel else config.OUTPUT_DIRECTORY
    directory.mkdir(parents=True, exist_ok=True)
    name = config.PANEL_NAMES[panel] if panel else config.FIGURE_NAME
    stem = directory / name
    document = pymupdf.open(stream=data, filetype="pdf")
    page = document[0]
    if panel:
        if panel == "b":
            page.draw_rect(pymupdf.Rect(350, 0, 367, 159), color=None, fill=(1, 1, 1), overlay=True)
        if panel == "c":
            page.draw_rect(pymupdf.Rect(359, 155, 379, 173), color=None, fill=(1, 1, 1), overlay=True)
            page.draw_rect(pymupdf.Rect(650, 45, 657, 192), color=None, fill=(1, 1, 1), overlay=True)
        if panel == "d":
            page.draw_rect(pymupdf.Rect(651, 12, 662, 43), color=None, fill=(1, 1, 1), overlay=True)
        if panel == "e":
            page.draw_rect(pymupdf.Rect(320, 172, 480, 200), color=None, fill=(1, 1, 1), overlay=True)
        page.set_cropbox(pymupdf.Rect(config.PANEL_BOXES[panel]))
        document.save(stem.with_suffix(".pdf"), garbage=3, deflate=True)
    else:
        stem.with_suffix(".pdf").write_bytes(data)
    stem.with_suffix(".svg").write_text(compact_svg(page.get_svg_image(text_as_path=True)), encoding="utf-8")
    page.get_pixmap(dpi=config.PNG_DPI, alpha=False).save(stem.with_suffix(".png"))
    document.close()
    print(f"Saved {stem.name}: PDF, SVG, PNG", flush=True)


def render_panel(panel):
    export_artifacts(assembled_pdf(), panel)


def render_figure():
    data = assembled_pdf()
    for panel in config.PANEL_NAMES:
        export_artifacts(data, panel)
    export_artifacts(data)
