"""Resolve SVG clipping to vector geometry and embedded image pixels."""

import base64
from copy import deepcopy
from functools import lru_cache
from io import BytesIO
import math
import re
import xml.etree.ElementTree as ET

from fontTools.misc.transform import Identity
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.svgLib.path import parse_path
import pathops
from PIL import Image, ImageChops

from .standalone_svg import (
    SVG, XLINK, URL, transform, mapped, path_bounds, intersection,
    style_properties, values, references, rectangular_path,
)

PAINT_KEYS = {
    "fill", "fill-opacity", "fill-rule", "stroke", "stroke-opacity", "stroke-width",
    "stroke-linecap", "stroke-linejoin", "stroke-miterlimit", "stroke-dasharray",
    "stroke-dashoffset", "color",
}


def matrix_text(matrix):
    return "matrix(" + " ".join(f"{value:.12g}" for value in matrix) + ")"


def contains(outer, inner, tolerance=1e-7):
    return inner is not None and outer[0] <= inner[0] + tolerance and outer[1] <= inner[1] + tolerance and outer[2] >= inner[2] - tolerance and outer[3] >= inner[3] - tolerance


def rectangle_path(bounds):
    x0, y0, x1, y1 = bounds
    path = pathops.Path()
    path.moveTo(x0, y0)
    path.lineTo(x1, y0)
    path.lineTo(x1, y1)
    path.lineTo(x0, y1)
    path.close()
    return path


@lru_cache(maxsize=4096)
def parsed_path(data):
    path = pathops.Path()
    parse_path(data, path.getPen())
    return path


def shape_path(element):
    tag = element.tag.removeprefix(SVG)
    if tag == "path":
        return pathops.Path(parsed_path(element.get("d", "")))
    if tag == "rect":
        x, y = float(element.get("x", 0)), float(element.get("y", 0))
        if float(element.get("rx", 0)) or float(element.get("ry", 0)):
            raise ValueError("Rounded rectangles require explicit path geometry.")
        return rectangle_path((x, y, x + float(element.get("width", 0)), y + float(element.get("height", 0))))
    if tag in ("circle", "ellipse"):
        x, y = float(element.get("cx", 0)), float(element.get("cy", 0))
        rx, ry = float(element.get("rx", element.get("r", 0))), float(element.get("ry", element.get("r", 0)))
        data = f"M{x-rx},{y}A{rx},{ry} 0 1 0 {x+rx},{y}A{rx},{ry} 0 1 0 {x-rx},{y}Z"
        return pathops.Path(parsed_path(data))
    path = pathops.Path()
    if tag == "line":
        path.moveTo(float(element.get("x1", 0)), float(element.get("y1", 0)))
        path.lineTo(float(element.get("x2", 0)), float(element.get("y2", 0)))
        return path
    if tag in ("polygon", "polyline"):
        points = values(element.get("points", ""))
        if points:
            path.moveTo(*points[:2])
            for index in range(2, len(points), 2):
                path.lineTo(*points[index:index + 2])
            if tag == "polygon":
                path.close()
        return path
    raise ValueError(f"Unsupported clipped shape: {tag}")


def image_data(element):
    href = element.get(XLINK, element.get("href", ""))
    if not href.startswith("data:image/png;base64,"):
        raise ValueError("Image must be an embedded PNG.")
    return Image.open(BytesIO(base64.b64decode(href.split(",", 1)[1]))).copy()


def set_image_data(element, image):
    output = BytesIO()
    image.save(output, format="PNG", optimize=True)
    element.attrib.pop("href", None)
    element.set(XLINK, "data:image/png;base64," + base64.b64encode(output.getvalue()).decode("ascii"))


def unmasked_svg(svg, return_report=False):
    ET.register_namespace("", SVG[1:-1])
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    root = ET.fromstring(svg)
    ids = {element.get("id"): element for element in root.iter() if element.get("id")}
    report = {"resolved_image_masks": 0, "removed_clipping_references": 0,
              "clipped_vector_elements": 0, "cropped_images": 0}
    for element in list(root.iter()):
        match = URL.fullmatch(element.get("mask", ""))
        if not match:
            continue
        mask = ids[match[1]]
        mask_images, source_images = list(mask.iter(SVG + "image")), list(element.iter(SVG + "image"))
        if len(mask_images) != 1 or len(source_images) != 1 or len(mask) != 1 or len(element) != 1:
            raise ValueError("Unsupported mask structure.")
        if transform(mask[0].get("transform", "")) != transform(element[0].get("transform", "")):
            raise ValueError("Mask and image transformations differ.")
        source = image_data(source_images[0]).convert("RGBA")
        alpha = image_data(mask_images[0]).convert("L")
        if source.size != alpha.size:
            raise ValueError("Mask and image dimensions differ.")
        source.putalpha(ImageChops.multiply(source.getchannel("A"), alpha))
        set_image_data(source_images[0], source)
        element.attrib.pop("mask")
        report["resolved_image_masks"] += 1

    defaults = {"fill": "black", "stroke": "none", "stroke-width": "1", "stroke-linecap": "butt", "stroke-linejoin": "miter"}
    for style in root.iter(SVG + "style"):
        if "stroke-linejoin: round" in (style.text or ""):
            defaults["stroke-linejoin"] = "round"

    def paint_for(element, inherited):
        return inherited | {key: value for key, value in style_properties(element).items() if key in PAINT_KEYS}

    def clip_box(element, matrix):
        match = URL.fullmatch(element.get("clip-path", ""))
        if not match:
            return None
        clip = ids[match[1]]
        if clip.get("clipPathUnits", "userSpaceOnUse") != "userSpaceOnUse" or len(clip) != 1:
            raise ValueError("Unsupported clipping units or compound clip.")
        child = clip[0]
        matrix = matrix.transform(transform(clip.get("transform", ""))).transform(transform(child.get("transform", "")))
        if abs(matrix[1]) + abs(matrix[2]) > 1e-9:
            raise ValueError("Clip axes must be aligned with the page.")
        if child.tag == SVG + "rect":
            x, y = float(child.get("x", 0)), float(child.get("y", 0))
            bounds = x, y, x + float(child.get("width", 0)), y + float(child.get("height", 0))
        elif child.tag == SVG + "path" and rectangular_path(child.get("d", "")):
            bounds = path_bounds(child.get("d"))
        else:
            raise ValueError("Clip must be a rectangle.")
        return mapped(bounds, matrix)

    def vector_parts(element, matrix, paint):
        if element.tag == SVG + "use":
            reference = element.get(XLINK, element.get("href", ""))
            target = ids[reference[1:]]
            matrix = matrix.translate(float(element.get("x", 0)), float(element.get("y", 0))).transform(transform(target.get("transform", "")))
            return vector_parts(target, matrix, paint_for(target, paint))
        if element.tag in (SVG + "g", SVG + "symbol"):
            parts = []
            for child in element:
                parts.extend(vector_parts(child, matrix.transform(transform(child.get("transform", ""))), paint_for(child, paint)))
            return parts
        path = shape_path(element)
        parts = []
        if paint.get("fill", "black") != "none":
            filled = pathops.Path(path)
            filled.fillType = pathops.FillType.EVEN_ODD if paint.get("fill-rule") == "evenodd" else pathops.FillType.WINDING
            parts.append((filled.transform(*matrix), paint.get("fill", "black"), paint.get("fill-opacity", "1")))
        if paint.get("stroke", "none") != "none":
            stroke = pathops.Path(path)
            cap = {"butt": pathops.LineCap.BUTT_CAP, "round": pathops.LineCap.ROUND_CAP, "square": pathops.LineCap.SQUARE_CAP}[paint.get("stroke-linecap", "butt")]
            join = {"miter": pathops.LineJoin.MITER_JOIN, "round": pathops.LineJoin.ROUND_JOIN, "bevel": pathops.LineJoin.BEVEL_JOIN}[paint.get("stroke-linejoin", "miter")]
            dash = values(paint.get("stroke-dasharray", "")) or None
            if dash and len(dash) % 2:
                dash *= 2
            stroke.stroke(float(paint.get("stroke-width", "1")), cap, join,
                          float(paint.get("stroke-miterlimit", "4")), dash, float(paint.get("stroke-dashoffset", 0)))
            stroke.convertConicsToQuads(0.0001)
            parts.append((stroke.transform(*matrix), paint["stroke"], paint.get("stroke-opacity", "1")))
        return parts

    def crop_image(element, matrix, parent_matrix, clip):
        source = image_data(element).convert("RGBA")
        pixel_matrix = matrix.translate(float(element.get("x", 0)), float(element.get("y", 0))).scale(float(element.get("width")) / source.width, float(element.get("height")) / source.height)
        if abs(pixel_matrix[1]) + abs(pixel_matrix[2]) > 1e-9:
            raise ValueError("Image axes must be aligned with the page.")
        extent = intersection(mapped(clip, pixel_matrix.inverse()), (0, 0, source.width, source.height))
        if extent is None:
            return None
        if contains(extent, (0, 0, source.width, source.height)):
            return element
        left, top, right, bottom = extent
        width, height = max(1, math.ceil(right-left)), max(1, math.ceil(bottom-top))
        cropped = source.transform((width, height), Image.Transform.EXTENT, extent, resample=Image.Resampling.BICUBIC)
        set_image_data(element, cropped)
        element.set("x", "0")
        element.set("y", "0")
        element.set("width", str(width))
        element.set("height", str(height))
        element.set("preserveAspectRatio", "none")
        local = parent_matrix.inverse().transform(pixel_matrix).translate(left, top).scale((right-left)/width, (bottom-top)/height)
        element.set("transform", matrix_text(local))
        report["cropped_images"] += 1
        return element

    def visit(element, parent_matrix, inherited, inherited_clip):
        tag = element.tag.removeprefix(SVG)
        if tag in ("defs", "metadata", "style", "title", "desc"):
            return element
        matrix = parent_matrix.transform(transform(element.get("transform", "")))
        paint = paint_for(element, inherited)
        own_clip = clip_box(element, matrix)
        clip = inherited_clip
        if own_clip is not None:
            clip = own_clip if clip is None else intersection(clip, own_clip)
            element.attrib.pop("clip-path")
            report["removed_clipping_references"] += 1
            if clip is None:
                return None
        if tag in ("g", "svg"):
            children = [visit(child, matrix, paint, clip) for child in element]
            element[:] = [child for child in children if child is not None]
            return element if len(element) or tag == "svg" else None
        if clip is None:
            return element
        if tag == "image":
            return crop_image(element, matrix, parent_matrix, clip)
        parts = vector_parts(element, matrix, paint)
        if all(contains(clip, path.bounds) for path, _, _ in parts):
            return element
        report["clipped_vector_elements"] += 1
        container = ET.Element(SVG + "g")
        source_properties = style_properties(element)
        for key in ("id", "opacity"):
            if source_properties.get(key) is not None:
                container.set(key, source_properties[key])
        clip_path = rectangle_path(clip)
        for path, fill, opacity in parts:
            cropped = pathops.op(path, clip_path, pathops.PathOp.INTERSECTION)
            if not len(cropped):
                continue
            cropped = cropped.transform(*parent_matrix.inverse())
            cropped.convertConicsToQuads(0.0001)
            pen = SVGPathPen(None, ntos=lambda value: f"{value:.9g}")
            cropped.draw(pen)
            attrs = {"d": pen.getCommands(), "fill": paint.get("color", "black") if fill == "currentColor" else fill, "stroke": "none", "fill-rule": "nonzero"}
            if opacity != "1":
                attrs["fill-opacity"] = opacity
            ET.SubElement(container, SVG + "path", attrs)
        return container if len(container) else None

    visit(root, Identity, defaults, None)
    for parent in list(root.iter()):
        parent[:] = [child for child in parent if child.tag not in (SVG + "mask", SVG + "clipPath")]
    used = set()
    for child in root:
        if child.tag != SVG + "defs":
            used.update(references(child))
    pending = list(used)
    while pending:
        identifier = pending.pop()
        if identifier in ids:
            for reference in references(ids[identifier]) - used:
                used.add(reference)
                pending.append(reference)
    for definitions in root.iter(SVG + "defs"):
        definitions[:] = [child for child in definitions if child.tag == SVG + "style" or {e.get("id") for e in child.iter()} & used]

    def flatten(element):
        children = []
        for child in element:
            if child.tag == SVG + "g":
                flatten(child)
                if not child.attrib:
                    children.extend(child)
                    continue
            if child.tag == SVG + "defs" and not len(child):
                continue
            children.append(child)
        element[:] = children
    flatten(root)
    forbidden = [element for element in root.iter() if element.tag in (SVG + "mask", SVG + "clipPath") or element.get("clip-path") or element.get("mask")]
    if forbidden:
        raise ValueError("Unresolved SVG clipping or masking remains.")
    report["remaining_masks"] = 0
    report["remaining_clipping_paths"] = 0
    result = ET.tostring(root, encoding="unicode", xml_declaration=True)
    return (result, report) if return_report else result
