"""Place the supplied scatter SVG components in the Figure 2 layout."""

from copy import deepcopy
import hashlib
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import pymupdf
import figure2_config as config

DIRECTORY = config.METHODS_DIRECTORY / "shared" / "assets" / "figure2" / "panel_e_source"
SOURCE = DIRECTORY / "score_top100_grid_3x2_fitted.svg"
NS = "{http://www.w3.org/2000/svg}"
XLINK = "{http://www.w3.org/1999/xlink}href"
TARGETS = (
    (27.1299, 204.0359, 154.0873, 275.6326),
    (27.4213, 307.1350, 154.3787, 377.1663),
    (185.7643, 203.4669, 312.7217, 273.4981),
    (187.1286, 306.0020, 314.0860, 376.0333),
    (343.7549, 202.5670, 470.7123, 272.5982),
    (344.8036, 306.1699, 471.7610, 376.2011),
)


def remove_previous_scatter(content):
    start, end = content.index(b"/Fm217 Do"), content.index(b"/Fm226 Do")
    block = content[start:end]
    pattern = rb"q 1 0 0 1 [-\d.]+ [-\d.]+ cm\n0 0 m\n(?:(?!\nQ).)*?\n[ fBS]+\nQ"
    signatures = {
        b"0.118 0 0.232 0.047 0.316 0.131 c",
        b"0.271 0 0.532 0.108 0.724 0.3 c",
        b"0.38 0 0.745 0.151 1.013 0.42 c",
        b"0.338 0 0.663 0.134 0.902 0.374 c",
    }
    removed = []

    def replace(match):
        if match[0].splitlines()[2] in signatures:
            removed.append(match[0])
            return b""
        return match[0]

    block = re.sub(pattern, replace, block, flags=re.S)
    if len(removed) != 653:
        raise ValueError(f"Expected 653 scatter and legend paths; found {len(removed)}.")
    for number in (217, 221, 222, 223, 224, 225):
        block = block.replace(f"/Fm{number} Do".encode(), b"")
    return content[:start] + block + content[end:]


def scatter_overlay():
    ET.register_namespace("", NS[1:-1])
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    source = ET.parse(SOURCE).getroot()
    root = ET.Element(NS + "svg", {"width": "841.89pt", "height": "595.276pt", "viewBox": "0 0 841.89 595.276"})
    definitions = ET.SubElement(root, NS + "defs")
    ids = set()
    for definition_group in source.iter(NS + "defs"):
        for item in definition_group:
            if item.get("id") not in ids:
                definitions.append(deepcopy(item))
                ids.add(item.get("id"))
    green = next(item for item in definitions if item.get("id") == "m3e090f9cd4")
    circle = deepcopy(green)
    circle.set("id", "delta6_open_circle")
    circle.attrib.pop("style", None)
    definitions.append(circle)
    counts = {"open_black": 0, "green": 0, "images": 0}
    for index, target in enumerate(TARGETS, 1):
        axes = next(item for item in source.iter(NS + "g") if item.get("id") == f"axes_{index}")
        coordinates = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", axes[0][0].get("d"))]
        xs, ys = coordinates[::2], coordinates[1::2]
        x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
        left, top, right, bottom = target
        sx, sy = (right - left) / (x1 - x0), (bottom - top) / (y1 - y0)
        clip_id = f"panel_e_plot_{index}"
        clip = ET.SubElement(definitions, NS + "clipPath", {"id": clip_id})
        ET.SubElement(clip, NS + "rect", {"x": str(left + 0.3), "y": str(top + 0.3), "width": str(right - left - 0.6), "height": str(bottom - top - 0.6)})
        outer = ET.SubElement(root, NS + "g", {"clip-path": f"url(#{clip_id})"})
        group = ET.SubElement(outer, NS + "g", {"transform": f"matrix({sx} 0 0 {sy} {left-sx*x0} {top-sy*y0})"})
        for original in axes:
            if original.tag != NS + "image" and not original.get("id", "").startswith("PathCollection_"):
                continue
            item = deepcopy(original)
            if item.tag == NS + "image":
                counts["images"] += 1
            for marker in item.iter(NS + "use"):
                style = marker.get("style", "")
                if "#d62728" in style:
                    marker.set(XLINK, "#delta6_open_circle")
                    marker.set("style", "fill: none; stroke: #000000; stroke-width: 0.8")
                    counts["open_black"] += 1
                elif "#00a359" in style:
                    marker.set("style", style.replace("#00a359", "#00a650"))
                    counts["green"] += 1
            group.append(item)
    ET.SubElement(root, NS + "circle", {"cx": "97.0375", "cy": "193.1996", "r": "1.64025", "fill": "none", "stroke": "#000000", "stroke-width": "0.4"})
    if counts != {"open_black": 592, "green": 30, "images": 6}:
        raise ValueError(f"Unexpected source SVG components: {counts}")
    svg = ET.tostring(root, encoding="utf-8", xml_declaration=True)
    svg_path = DIRECTORY / "panel_e_open_circle_components.svg"
    pdf_path = svg_path.with_suffix(".pdf")
    hash_path = svg_path.with_suffix(".sha256")
    digest = hashlib.sha256(svg).hexdigest()
    svg_path.write_bytes(svg)
    if not pdf_path.exists() or not hash_path.exists() or hash_path.read_text().strip() != digest:
        converter = shutil.which("rsvg-convert") or "/opt/anaconda3/bin/rsvg-convert"
        if not Path(converter).is_file():
            raise FileNotFoundError("Install librsvg to render a changed source SVG.")
        subprocess.run([converter, "--format=pdf", "--output", str(pdf_path), str(svg_path)], check=True)
        hash_path.write_text(digest + "\n")
    return pdf_path


def add_scatter_components(data):
    document = pymupdf.open(stream=data, filetype="pdf")
    with pymupdf.open(scatter_overlay()) as overlay:
        document[0].show_pdf_page(document[0].rect, overlay, 0)
    result = document.tobytes(garbage=3, deflate=True)
    document.close()
    return result
