"""Remove SVG artwork outside its page and discard unused definitions."""

from functools import lru_cache
import math
import re
import xml.etree.ElementTree as ET

from fontTools.misc.transform import Transform, Identity
from fontTools.pens.boundsPen import BoundsPen
from fontTools.pens.recordingPen import RecordingPen
from fontTools.svgLib.path import parse_path

SVG = "{http://www.w3.org/2000/svg}"
XLINK = "{http://www.w3.org/1999/xlink}href"
NUMBERS = re.compile(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?")
URL = re.compile(r"url\(\s*['\"]?#([^)'\"\s]+)['\"]?\s*\)")


def values(text):
    return [float(item) for item in NUMBERS.findall(text or "")]


@lru_cache(maxsize=65536)
def transform(text):
    result = Identity
    for name, arguments in re.findall(r"([A-Za-z]+)\s*\(([^)]*)\)", text or ""):
        args = values(arguments)
        if name == "matrix":
            item = Transform(*args)
        elif name == "translate":
            item = Identity.translate(args[0], args[1] if len(args) > 1 else 0)
        elif name == "scale":
            item = Identity.scale(args[0], args[1] if len(args) > 1 else args[0])
        elif name == "rotate":
            angle = math.radians(args[0])
            item = Identity.rotate(angle)
            if len(args) == 3:
                item = Identity.translate(args[1], args[2]).transform(item).translate(-args[1], -args[2])
        elif name == "skewX":
            item = Transform(1, 0, math.tan(math.radians(args[0])), 1, 0, 0)
        elif name == "skewY":
            item = Transform(1, math.tan(math.radians(args[0])), 0, 1, 0, 0)
        else:
            raise ValueError(f"Unsupported SVG transform: {name}")
        result = result.transform(item)
    return result


@lru_cache(maxsize=65536)
def path_bounds(path):
    pen = BoundsPen(None)
    parse_path(path, pen)
    return pen.bounds


@lru_cache(maxsize=65536)
def rectangular_path(path):
    pen = RecordingPen()
    parse_path(path, pen)
    if any(command not in ("moveTo", "lineTo", "closePath", "endPath") for command, _ in pen.value):
        return False
    points = [point for command, args in pen.value for point in args]
    if len(points) not in (4, 5):
        return False
    xx, yy = {point[0] for point in points}, {point[1] for point in points}
    return len(xx) == len(yy) == 2 and set(points) == {(x, y) for x in xx for y in yy}


def union(bounds):
    bounds = [item for item in bounds if item is not None]
    if not bounds:
        return None
    return min(item[0] for item in bounds), min(item[1] for item in bounds), max(item[2] for item in bounds), max(item[3] for item in bounds)


def intersection(first, second):
    if first is None or second is None:
        return None
    result = max(first[0], second[0]), max(first[1], second[1]), min(first[2], second[2]), min(first[3], second[3])
    return result if result[0] <= result[2] and result[1] <= result[3] else None


def contains(outer, inner, tolerance=0.00001):
    return inner is not None and outer[0] - tolerance <= inner[0] and outer[1] - tolerance <= inner[1] and outer[2] + tolerance >= inner[2] and outer[3] + tolerance >= inner[3]


def mapped(bounds, matrix):
    if bounds is None:
        return None
    points = [matrix.transformPoint(point) for point in ((bounds[0], bounds[1]), (bounds[0], bounds[3]), (bounds[2], bounds[1]), (bounds[2], bounds[3]))]
    return min(point[0] for point in points), min(point[1] for point in points), max(point[0] for point in points), max(point[1] for point in points)


def style_properties(element):
    result = dict(element.attrib)
    for declaration in element.get("style", "").split(";"):
        if ":" in declaration:
            key, value = declaration.split(":", 1)
            result[key.strip()] = value.strip()
    return result


def references(element):
    result = set()
    for node in element.iter():
        for value in node.attrib.values():
            result.update(URL.findall(value))
        link = node.get(XLINK, node.get("href", ""))
        if link.startswith("#"):
            result.add(link[1:])
    return result


def standalone_svg(svg, return_report=False):
    ET.register_namespace("", SVG[1:-1])
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    root = ET.fromstring(svg)
    x, y, width, height = values(root.get("viewBox"))
    viewport = (x, y, x + width, y + height)
    ids = {element.get("id"): element for element in root.iter() if element.get("id")}
    report = {"original_elements": sum(1 for _ in root.iter()), "removed_graphics": 0, "removed_groups": 0, "removed_definitions": 0, "removed_covered_graphics": 0, "removed_redundant_clips": 0, "remaining_boundary_objects": []}
    local_cache = {}
    graphics = []
    kept_bounds = {}

    def local_bounds(element):
        key = id(element)
        if key in local_cache:
            return local_cache[key]
        tag = element.tag.removeprefix(SVG)
        props = style_properties(element)
        if tag == "path":
            bounds = path_bounds(element.get("d", ""))
        elif tag in ("rect", "image"):
            xx, yy = float(element.get("x", 0)), float(element.get("y", 0))
            bounds = (xx, yy, xx + float(element.get("width", 0)), yy + float(element.get("height", 0)))
        elif tag in ("circle", "ellipse"):
            xx, yy = float(element.get("cx", 0)), float(element.get("cy", 0))
            rx = float(element.get("rx", element.get("r", 0)))
            ry = float(element.get("ry", element.get("r", 0)))
            bounds = (xx - rx, yy - ry, xx + rx, yy + ry)
        elif tag == "line":
            xx = [float(element.get(name, 0)) for name in ("x1", "x2")]
            yy = [float(element.get(name, 0)) for name in ("y1", "y2")]
            bounds = (min(xx), min(yy), max(xx), max(yy))
        elif tag in ("polygon", "polyline"):
            points = values(element.get("points"))
            bounds = (min(points[::2]), min(points[1::2]), max(points[::2]), max(points[1::2])) if points else None
        elif tag == "use":
            link = element.get(XLINK, element.get("href", ""))
            if not link.startswith("#") or link[1:] not in ids:
                raise ValueError(f"Unresolved SVG reference: {link}")
            target = ids[link[1:]]
            bounds = mapped(local_bounds(target), Identity.translate(float(element.get("x", 0)), float(element.get("y", 0))).transform(transform(target.get("transform", ""))))
        elif tag in ("g", "clipPath", "mask", "symbol", "svg"):
            bounds = union(mapped(local_bounds(child), transform(child.get("transform", ""))) for child in element if child.tag != SVG + "defs")
        elif tag in ("defs", "metadata", "title", "desc", "style"):
            bounds = None
        else:
            raise ValueError(f"Unsupported SVG element: {tag}")
        if bounds is not None and props.get("stroke", "none") != "none":
            padding = float(props.get("stroke-width", "1")) / 2
            bounds = (bounds[0] - padding, bounds[1] - padding, bounds[2] + padding, bounds[3] + padding)
        local_cache[key] = bounds
        return bounds

    def clip_bounds(element, matrix):
        clip = element.get("clip-path", "")
        match = URL.fullmatch(clip)
        if match and match[1] in ids:
            definition = ids[match[1]]
            if definition.get("clipPathUnits", "userSpaceOnUse") == "userSpaceOnUse":
                return mapped(local_bounds(definition), matrix.transform(transform(definition.get("transform", ""))))
        return None

    def is_rectangle(element):
        if element.tag == SVG + "rect":
            return float(element.get("rx", 0)) == 0 and float(element.get("ry", 0)) == 0
        if element.tag == SVG + "path":
            return rectangular_path(element.get("d", ""))
        return len(element) == 1 and is_rectangle(element[0])

    def visit(element, parent_matrix, visible):
        tag = element.tag.removeprefix(SVG)
        if tag in ("defs", "metadata", "title", "desc", "style"):
            return True
        matrix = parent_matrix.transform(transform(element.get("transform", "")))
        props = style_properties(element)
        if props.get("display") == "none" or props.get("visibility") == "hidden" or props.get("opacity") == "0":
            report["removed_graphics"] += 1
            return False
        own_clip = clip_bounds(element, matrix)
        if own_clip is not None:
            visible = intersection(visible, own_clip)
        if visible is None:
            report["removed_groups"] += 1
            return False
        if tag in ("g", "svg"):
            children = [child for child in element if visit(child, matrix, visible)]
            element[:] = children
            if not children:
                report["removed_groups"] += 1
                return False
            kept_bounds[id(element)] = union(kept_bounds.get(id(child)) for child in children)
            match = URL.fullmatch(element.get("clip-path", ""))
            if match and own_clip is not None and is_rectangle(ids[match[1]]) and (contains(own_clip, viewport) or contains(own_clip, kept_bounds[id(element)])):
                element.attrib.pop("clip-path")
                report["removed_redundant_clips"] += 1
            return True
        bounds = mapped(local_bounds(element), matrix)
        if bounds is None or intersection(bounds, visible) is None:
            report["removed_graphics"] += 1
            return False
        if not contains(viewport, bounds) and is_rectangle(element) and props.get("fill", "").lower() in ("white", "#ffffff", "#fff") and abs(matrix[1]) + abs(matrix[2]) < 1e-12:
            clipped = mapped(intersection(bounds, viewport), matrix.inverse())
            element.tag = SVG + "rect"
            element.attrib.pop("d", None)
            for key, value in zip(("x", "y", "width", "height"), (clipped[0], clipped[1], clipped[2]-clipped[0], clipped[3]-clipped[1])):
                element.set(key, f"{value:.12g}")
            bounds = intersection(bounds, viewport)
        kept_bounds[id(element)] = bounds
        graphics.append((element, bounds, intersection(bounds, visible), props, is_rectangle(element)))
        return True

    root[:] = [child for child in root if visit(child, Identity, viewport)]
    occluders, covered = [], set()
    for element, bounds, visible_bounds, props, rectangle in reversed(graphics):
        if any(contains(box, visible_bounds) for box in occluders):
            covered.add(id(element))
            report["removed_covered_graphics"] += 1
            continue
        if rectangle and props.get("fill", "").lower() in ("white", "#ffffff", "#fff") and float(props.get("fill-opacity", 1)) == 1 and float(props.get("opacity", 1)) == 1:
            occluders.append(visible_bounds)
    def remove_covered(element):
        kept = []
        for child in element:
            if id(child) in covered:
                continue
            if child.tag == SVG + "g":
                remove_covered(child)
                if not len(child):
                    continue
            kept.append(child)
        element[:] = kept
    remove_covered(root)
    for element, bounds, _, props, _ in graphics:
        if id(element) not in covered and not contains(viewport, bounds):
            report["remaining_boundary_objects"].append({"tag": element.tag.removeprefix(SVG), "bounds": list(bounds), "id": element.get("id"), "text": element.get("data-text"), "fill": props.get("fill"), "d": element.get("d", "")[:160]})
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
    for definitions in list(root.iter(SVG + "defs")):
        kept = []
        for element in definitions:
            element_ids = {node.get("id") for node in element.iter()}
            if element.tag == SVG + "style" or element_ids & used:
                kept.append(element)
            else:
                report["removed_definitions"] += 1
        definitions[:] = kept
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
    report["final_elements"] = sum(1 for _ in root.iter())
    result = ET.tostring(root, encoding="unicode", xml_declaration=True)
    return (result, report) if return_report else result
