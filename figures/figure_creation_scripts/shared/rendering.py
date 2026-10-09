"""Vector plotting and PDF, SVG, and PNG export for Figure 1."""

from io import BytesIO
import re
import xml.etree.ElementTree as ET

import numpy as np
import pymupdf
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject

from shared import analysis
import config

C_BASE = 377.2461093920554
C_SCALE = (431.2989661547525 - C_BASE) / 124
D_BASE = 397.7699429776292
D_SCALE = 0.661994
D_SET_BASE = 551.477167735118
D_SET_SCALE = 0.1394
INTERSECTION_ORDER = ((0, 0, 0, 1), (0, 0, 1, 1), (1, 0, 1, 1), (1, 1, 1, 1), (0, 1, 0, 0))


def numbers(values):
    return " ".join(f"{float(value):.12f}".rstrip("0").rstrip(".") if value else "0" for value in values)


def update_rectangles(content, counts, intersections):
    replacements = []
    for x, strain in zip((431.688653380506, 450.498987831803, 469.309317146986, 488.119651598281), config.STRAINS):
        replacements.append((x, C_BASE, 11.66240725708, counts[strain] * C_SCALE))
    for x, pattern in zip((592.452583870936, 604.18027508334, 615.907966295747, 627.63565750815, 639.363348720557), INTERSECTION_ORDER):
        replacements.append((x, D_BASE, 5.863845606202, intersections.get(pattern, 0) * D_SCALE))
    for y, strain in zip((351.754726561239, 362.12518617026, 372.495645779281, 382.86610903512), reversed(config.STRAINS)):
        replacements.append((D_SET_BASE, y, -counts[strain] * D_SET_SCALE, 5.185229804516))
    pattern = re.compile(rb"(?m)^(-?[\d.]+) (-?[\d.]+) (-?[\d.]+) (-?[\d.]+) re\n([BfS])")

    def replace(match):
        x, y = map(float, match.group(1, 2))
        for xx, yy, width, height in replacements:
            if abs(x - xx) < 0.00001 and abs(y - yy) < 0.00001:
                return numbers((xx, yy, width, height)).encode() + b" re\n" + match[5]
        return match[0]

    return pattern.sub(replace, content)


def text_command(size, x, y, text):
    return numbers((size, 0, 0, size, x, y)).encode() + b" Tm\n(" + str(text).encode("ascii") + b")Tj\n"


def replace_text(content, size, old_x, old_y, x, y, text):
    pattern = re.compile(rb"([\d.]+) 0 0 ([\d.]+) ([\d.]+) ([\d.]+) Tm\n\([^\n]*?\)Tj\n")
    hits = 0

    def replace(match):
        nonlocal hits
        sx, sy, xx, yy = map(float, match.group(1, 2, 3, 4))
        if max(abs(sx - size), abs(sy - size), abs(xx - old_x), abs(yy - old_y)) < 0.00001:
            hits += 1
            return text_command(size, x, y, text)
        return match[0]

    result = pattern.sub(replace, content)
    if hits != 1:
        raise ValueError(f"Expected one text anchor at ({old_x}, {old_y}); found {hits}.")
    return result


def update_hit_labels(content, counts, intersections):
    c_labels = (
        ("WT", 8, 435.154935, 382.480530135183),
        ("Delta40", 5, 453.60903, 380.730210135183),
        ("Delta3", 39, 470.699895, 395.427809729634),
        ("Delta6", 124, 486.654735, 433.485515),
    )
    for strain, previous, x, y in c_labels:
        value = counts[strain]
        xx = x + (len(str(previous)) - len(str(value))) * 0.556 * 8.415 / 2
        yy = y + (value - previous) * C_SCALE
        content = replace_text(content, 8.415, x, y, xx, yy, value)
    for strain, x, y in (("Delta3", 450.3945, 436.1699), ("Delta6", 455.788515, 448.85972)):
        value = f"~{round(counts[strain] / counts['WT'])}x" if counts["WT"] else "NA"
        content = replace_text(content, 8.415, x, y, x, y, value)
    d_top = (
        (592.2441, 456.149432977629),
        (604.12776, 420.401756977629),
        (616.746795, 402.527918977629),
        (628.965786, 402.527918977629),
        (640.514115, 400.541936977629),
    )
    for pattern, (x, y) in zip(INTERSECTION_ORDER, d_top):
        value = intersections.get(pattern, 0)
        content = replace_text(content, 5.883, x, y, x, D_BASE + value * D_SCALE + 2.11, value)
    for strain, x, y in (
        ("WT", 545.242807735118, 383.5229),
        ("Delta40", 545.661007735118, 373.4888),
        ("Delta3", 537.802247735118, 362.70638),
        ("Delta6", 522.834087735118, 352.2269),
    ):
        value = counts[strain]
        xx = D_SET_BASE - value * D_SET_SCALE - 2.0 - len(str(value)) * 0.556 * 5.61
        content = replace_text(content, 5.61, x, y, xx, y, value)
    start = content.index(b"8.211 0 0 8.211 801.4209 394.6567 Tm")
    end = content.index(b"0 0.418188750744 1 0  scn", start)
    x, y = 801.4209, 394.6567
    blocks = [text_command(8.211, x, y, intersections.get((0, 0, 0, 1), 0))]
    for dx, dy, pattern in (
        (-1.266, 3.234, (0, 0, 1, 1)),
        (-4.162, -3.442, (1, 1, 1, 1)),
        (-2.139, -2.331, (1, 0, 1, 1)),
    ):
        x += dx * 8.211
        y += dy * 8.211
        blocks.append(text_command(8.211, x, y, intersections.get(pattern, 0)))
    blocks.append(text_command(9.384, 728.5381, 439.7817, intersections.get((0, 1, 0, 0), 0)))
    return content[:start] + b"".join(blocks) + content[end:]


def updated_layout():
    page = PdfReader(config.ASSET_DIRECTORY / "figure_1_vector_layout.pdf").pages[0]
    counts, intersections = analysis.hit_summary()
    active_patterns = {pattern for pattern, count in intersections.items() if any(pattern) and count}
    if not active_patterns.issubset(INTERSECTION_ORDER):
        raise ValueError("Update the intersection layout for the supplied strain combinations.")
    content = page.get_contents().get_data()
    content = update_rectangles(content, counts, intersections)
    return page, update_hit_labels(content, counts, intersections)


def chemical_space_stream(panel):
    frame = analysis.panel_coordinates(panel)
    coordinates = frame[["tSNE_1", "tSNE_2"]].to_numpy()
    if not np.isfinite(coordinates).all():
        raise ValueError("Chemical-space coordinates must be finite.")
    bounds = ((366, 606), (87, 310)) if panel == "f" else ((631, 856), (91, 311))
    ranges = np.ptp(coordinates, axis=0)
    center = (coordinates.min(axis=0) + coordinates.max(axis=0)) / 2
    scale = min((bounds[index][1] - bounds[index][0]) / ranges[index] for index in range(2))
    xy = (coordinates - center) * scale + np.array([sum(value) / 2 for value in bounds])
    wt_columns = ["WT_8K", "WT_39K", "WT_ECBD"] if panel == "f" else ["WT_8K", "WT_COADD"]
    wt = frame[wt_columns].eq(1).any(axis=1).to_numpy()
    delta6 = frame.D6_8K.eq(1).to_numpy()
    delta5 = frame.D5_COADD.eq(1).to_numpy() if panel == "g" else np.zeros(len(frame), dtype=bool)
    background = ~(wt | delta6 | delta5)
    gray = (0.421875, 0.34375, 0.34765625, 0.0078125) if panel == "f" else (0.23046875, 0.1796875, 0.18359375, 0)
    radius = 1.691855768689 if panel == "f" else 1.705688905716
    layers = (
        ("background", background, gray, 0.30065052797 if panel == "f" else 0.303108739853),
        ("Delta5", delta5, config.COLORS_CMYK["Delta5"], 1.705688905716),
        ("Delta6", delta6, config.COLORS_CMYK["Delta6"], radius),
        ("WT", wt, config.COLORS_CMYK["WT"], radius),
    )
    commands = ["q", "/CS0 cs /CS0 CS", "0 J 0 j 0.25 w 4 M", "[] 0 d"]
    for name, selected, color, r in layers:
        if not selected.any():
            continue
        graphics_state = ("/GS6" if panel == "f" else "/GS3") if name == "background" else ("/GS0" if name == "WT" else "/GS2")
        commands.extend((numbers(color) + " scn", "0.74609375 0.67578125 0.66796875 0.8984375 SCN", graphics_state + " gs", f"{0.25 / r:.12f} w"))
        marker = "Fill" if name == "background" else "Stroke"
        for x, y in xy[selected]:
            commands.append(f"q {r:.12f} 0 0 {r:.12f} {x:.9f} {y:.9f} cm /FigDot{marker} Do Q")
    commands.append("Q")
    return ("\n".join(commands) + "\n").encode()


def pdf_bytes(page, content):
    stream = DecodedStreamObject()
    stream.set_data(content)
    page[NameObject("/Contents")] = stream.flate_encode()
    writer = PdfWriter()
    writer.add_page(page)
    writer.add_metadata({"/Title": "Antipseudomonal screening and chemical space"})
    buffer = BytesIO()
    writer.write(buffer)
    return buffer.getvalue()


def compact_svg(svg):
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    root = ET.fromstring(svg)
    ns = "{http://www.w3.org/2000/svg}"
    marker_clips = {}
    for element in root.iter(ns + "clipPath"):
        children = list(element)
        if len(children) == 1 and children[0].get("d") == "M-1.2-1.2H1.2V1.2H-1.2Z":
            marker_clips[element.get("id")] = children[0].get("transform")
    removed_clips = set()
    for group in root.iter(ns + "g"):
        reference = group.get("clip-path", "")
        identifier = reference[5:-1] if reference.startswith("url(#") else ""
        if identifier not in marker_clips:
            continue
        paths = list(group.iter(ns + "path"))
        if paths and all(
            path.get("d", "").startswith("M0-1C.2652165")
            and path.get("transform") == marker_clips[identifier]
            and float(path.get("stroke-width", "0")) <= 0.4
            for path in paths
        ):
            del group.attrib["clip-path"]
            removed_clips.add(identifier)
    for definitions_element in root.iter(ns + "defs"):
        definitions_element[:] = [child for child in definitions_element if child.get("id") not in removed_clips]
    definitions = ET.Element(ns + "defs")
    shapes = {}
    for element in list(root.iter(ns + "path")):
        path = element.get("d")
        if path is None or len(path) < 80 or element.get("id"):
            continue
        if path not in shapes:
            name = f"shared_marker_{len(shapes)}"
            shapes[path] = name
            ET.SubElement(definitions, ns + "path", {"id": name, "d": path})
        element.tag = ns + "use"
        del element.attrib["d"]
        element.set("{http://www.w3.org/1999/xlink}href", "#" + shapes[path])
    _, _, width, height = root.get("viewBox").split()
    canvas_clip = ET.SubElement(definitions, ns + "clipPath", {"id": "figure_page_clip"})
    ET.SubElement(canvas_clip, ns + "rect", {"x": "0", "y": "0", "width": width, "height": height})
    artwork = ET.Element(ns + "g", {"clip-path": "url(#figure_page_clip)"})
    existing_definitions = []
    for child in root:
        if child.tag == ns + "defs":
            existing_definitions.append(child)
        else:
            artwork.append(child)
    root[:] = existing_definitions
    root.insert(0, definitions)
    ET.SubElement(root, ns + "rect", {"x": "0", "y": "0", "width": width, "height": height, "fill": "#ffffff"})
    root.append(artwork)
    from shared.standalone_svg import standalone_svg
    from shared.unmasked_svg import unmasked_svg
    return unmasked_svg(standalone_svg(ET.tostring(root, encoding="unicode", xml_declaration=True)))


def export_artifacts(data, name, panel=None):
    output_directory = config.PANEL_DIRECTORIES[panel] if panel else config.OUTPUT_DIRECTORY
    output_directory.mkdir(parents=True, exist_ok=True)
    stem = output_directory / name
    document = pymupdf.open(stream=data, filetype="pdf")
    page = document[0]
    if panel:
        x0, y0, x1, y1 = config.PANEL_BOXES[panel]
        height = page.mediabox.height
        clip = pymupdf.Rect(x0, max(0, height - y1), x1, height - y0)
        page.set_cropbox(clip)
        document.save(stem.with_suffix(".pdf"), garbage=3, deflate=True)
    else:
        stem.with_suffix(".pdf").write_bytes(data)
    svg = page.get_svg_image(text_as_path=True)
    stem.with_suffix(".svg").write_text(compact_svg(svg), encoding="utf-8")
    page.get_pixmap(dpi=config.PNG_DPI, alpha=False).save(stem.with_suffix(".png"))
    document.close()
    print(f"Saved {stem.name}: PDF, SVG, PNG", flush=True)


def render_panel(panel):
    if panel in ("a", "b"):
        page = PdfReader(config.ASSET_DIRECTORY / "figure_1_vector_layout.pdf").pages[0]
        content = page.get_contents().get_data()
    else:
        page, content = updated_layout()
    if panel != "b":
        content = re.sub(rb"/Fm(\d+) Do", lambda match: b"" if int(match[1]) >= 27 else match[0], content)
        forms = page["/Resources"]["/XObject"]
        for name in list(forms):
            if name.startswith("/Fm") and int(name[3:]) >= 27:
                del forms[name]
    if panel in ("f", "g"):
        content += chemical_space_stream(panel)
    if panel == "b":
        from shared.figure1_insets import inset_stream
        content += inset_stream(page)
    export_artifacts(pdf_bytes(page, content), config.PANEL_NAMES[panel], panel)
    if panel not in ("a", "b"):
        analysis.write_summaries()


def render_figure():
    page, content = updated_layout()
    content += chemical_space_stream("f") + chemical_space_stream("g")
    from shared.figure1_insets import inset_stream
    content += inset_stream(page)
    export_artifacts(pdf_bytes(page, content), config.FIGURE_NAME)
    analysis.write_summaries()
