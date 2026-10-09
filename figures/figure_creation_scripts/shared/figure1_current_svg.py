"""Standalone SVG export preserving complex native density curves."""
import xml.etree.ElementTree as ET
from pathops import PathOpsError

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
    cleaned = standalone_svg(ET.tostring(root, encoding="unicode", xml_declaration=True))
    try:
        return unmasked_svg(cleaned)
    except PathOpsError:
        # The saved density artwork has complex paths that cannot always be
        # boolean-clipped by Skia. Standard SVG clipping preserves their exact
        # appearance. Off-page/covered objects have already been discarded.
        return cleaned

