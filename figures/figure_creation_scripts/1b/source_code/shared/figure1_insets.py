"""Vector zoom insets for the four replicate scatter plots."""

from functools import lru_cache
import json
import numpy as np
import pandas as pd
import pymupdf
from pypdf.generic import DictionaryObject, FloatObject, NameObject

import config

AXES = {
    "WT": (41.0019, 137.2154, 174.9550, 271.1684),
    "Delta40": (211.2520, 137.2154, 345.2050, 271.1684),
    "Delta3": (41.0019, 306.5778, 174.9550, 440.5308),
    "Delta6": (211.2520, 306.5778, 345.2050, 440.5308),
}
RGB = {
    "WT": (0.0, 0.649897, 0.315145),
    "Delta40": (0.979049, 0.637232, 0.102022),
    "Delta3": (0.752651, 0.228199, 0.169512),
    "Delta6": (0.658030, 0.377462, 0.647150),
}
INSET_MIN = 75.0
INSET_MAX = 105.0
INSET_SIZE = 44.0
POINT_RADIUS = 0.82
MINIMUM_CENTER_DISTANCE = 2.0
ARTWORK_OFFSET_X = -1.0
ARTWORK_OFFSET_Y = -0.726


def separate_points(coordinates):
    original = np.asarray(coordinates, dtype=float)
    if len(original) < 2:
        return original.copy()
    positions = original.copy()
    n = len(positions)
    first, second = np.triu_indices(n, 1)
    margin = POINT_RADIUS + 0.06
    for iteration in range(4000):
        displacement = positions[first] - positions[second]
        distance = np.linalg.norm(displacement, axis=1)
        close = distance < MINIMUM_CENTER_DISTANCE + 0.0001
        if not close.any():
            break
        i, j = first[close], second[close]
        vector = displacement[close]
        length = distance[close]
        zero = length < 1e-10
        if zero.any():
            angle = (i[zero] * 2.399963229728653 + j[zero] * 0.754877666246693) % (2 * np.pi)
            vector[zero] = np.column_stack((np.cos(angle), np.sin(angle)))
            length[zero] = 1
        push = vector / length[:, None] * ((MINIMUM_CENTER_DISTANCE + 0.0002 - distance[close]) / 2)[:, None]
        total = np.zeros_like(positions)
        np.add.at(total, i, push)
        np.add.at(total, j, -push)
        degree = np.bincount(np.r_[i, j], minlength=n)
        total /= np.maximum(1, degree)[:, None]
        positions += total * 0.9
        if iteration < 1200:
            positions += (original - positions) * 0.002
        positions = np.clip(positions, margin, INSET_SIZE - margin)
    distances = np.linalg.norm(positions[first] - positions[second], axis=1)
    if distances.min() < MINIMUM_CENTER_DISTANCE - 0.00001:
        raise ValueError("Inset marker spacing did not converge.")
    return positions


@lru_cache(maxsize=1)
def plotted_inset_points():
    source = pymupdf.open(config.ASSET_DIRECTORY / "figure_1_vector_layout.pdf")
    drawings = source[0].get_drawings()
    rows = []
    for strain, (x0, y0, x1, y1) in AXES.items():
        scale = (x1 - x0) / 166
        for item in drawings:
            if item["fill"] is None or not np.allclose(item["fill"], RGB[strain], rtol=0, atol=0.00001):
                continue
            rect = item["rect"]
            if not 1.3 < rect.width < 1.4 or not 1.3 < rect.height < 1.4:
                continue
            x, y = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
            if not (x0 <= x <= x1 and y0 <= y <= y1):
                continue
            first = (x - x0) / scale - 63
            second = (y1 - y) / scale - 63
            if INSET_MIN <= first <= INSET_MAX and INSET_MIN <= second <= INSET_MAX:
                rows.append({"Strain": strain, "Replicate_1_Growth_Inhibition": first, "Replicate_2_Growth_Inhibition": second})
    source.close()
    frame = pd.DataFrame(rows)
    report = {}
    for strain in AXES:
        selected = frame.Strain.eq(strain)
        values = frame.loc[selected, ["Replicate_1_Growth_Inhibition", "Replicate_2_Growth_Inhibition"]].to_numpy()
        original = (values - INSET_MIN) * INSET_SIZE / (INSET_MAX - INSET_MIN)
        display = separate_points(original) if strain == "WT" else original.copy()
        displayed_values = display * (INSET_MAX - INSET_MIN) / INSET_SIZE + INSET_MIN if strain == "WT" else values.copy()
        frame.loc[selected, "Display_Replicate_1_Growth_Inhibition"] = displayed_values[:, 0]
        frame.loc[selected, "Display_Replicate_2_Growth_Inhibition"] = displayed_values[:, 1]
        shift = np.linalg.norm(display - original, axis=1)
        frame.loc[selected, "Display_Offset_Distance_Points"] = shift
        distance = np.linalg.norm(display[:, None] - display[None, :], axis=2)
        np.fill_diagonal(distance, np.inf)
        report[strain] = {"points": len(values), "display_separation_applied": strain == "WT", "maximum_display_offset_points": float(shift.max()), "minimum_marker_center_distance_points": float(distance.min()), "marker_diameter_points": 2 * POINT_RADIUS}
    directory = config.PANEL_DIRECTORIES["b"]
    frame.to_csv(directory / "replicate_scatter_75_to_105_inset_coordinates.csv", index=False)
    (directory / "inset_display_validation.json").write_text(json.dumps({"axis_limits": [INSET_MIN, INSET_MAX], "strains": report, "main_plot_coordinates_changed": False}, indent=2) + "\n")
    return frame


def inset_stream(page):
    points = plotted_inset_points()
    height = float(page.mediabox.height)
    states = page["/Resources"]["/ExtGState"]
    states[NameObject("/InsetPointAlpha")] = DictionaryObject({NameObject("/ca"): FloatObject(0.7), NameObject("/CA"): FloatObject(1)})
    commands = [f"q 1 0 0 1 {ARTWORK_OFFSET_X} {-ARTWORK_OFFSET_Y} cm /CS0 cs /CS0 CS /GS0 gs 0 J 0 j [] 0 d"]

    def line(xa, ya, xb, yb):
        commands.append(f"{xa:.6f} {height-ya:.6f} m {xb:.6f} {height-yb:.6f} l S")

    def rectangle(x, y, w, h, operation):
        commands.append(f"{x:.6f} {height-y-h:.6f} {w:.6f} {h:.6f} re {operation}")

    def text(x, y, value, align="center", size=4.5):
        width = pymupdf.get_text_length(str(value), fontname="helv", fontsize=size)
        start = x - (width / 2 if align == "center" else width)
        commands.append(f"BT /UpdatedTT2 1 Tf 0 Tc 0 Tw 0 Ts 100 Tz 0 Tr {size} 0 0 {size} {start:.6f} {height-y:.6f} Tm ({value}) Tj ET")

    for strain, (x0, y0, x1, y1) in AXES.items():
        scale = (x1 - x0) / 166
        roi_x = x0 + (INSET_MIN + 63) * scale
        roi_y = max(y0, y1 - (INSET_MAX + 63) * scale)
        roi_right = min(x1, x0 + (INSET_MAX + 63) * scale)
        roi_bottom = y1 - (INSET_MIN + 63) * scale
        left, top, size = x1 - 49, y1 - 55, INSET_SIZE
        commands.extend(("0 0 0 0.48 SCN 0.3 w [1.6 1.6] 0 d",))
        rectangle(roi_x, roi_y, roi_right - roi_x, roi_bottom - roi_y, "S")
        line(roi_x, roi_bottom, left, top)
        line(roi_right, roi_bottom, left + size, top)
        commands.extend(("0 0 0 0 scn [] 0 d",))
        rectangle(left - 10.5, top - 9, size + 14, size + 18, "f")
        if strain == "WT":
            commands.append("0 0 0 0.65 scn")
            text(left + size / 2, top - 3, "points offset", size=4)
        commands.extend(("q",))
        rectangle(left, top, size, size, "W n")
        commands.extend(("0 0 0 0.65 SCN 0.4 w [2.2 1.8] 0 d",))
        line(left, top + size, left + size, top)
        commands.extend(("[] 0 d /InsetPointAlpha gs", " ".join(str(v) for v in config.COLORS_CMYK[strain]) + " scn"))
        for row in points.loc[points.Strain.eq(strain)].itertuples(index=False):
            x = left + (row.Display_Replicate_1_Growth_Inhibition - INSET_MIN) * size / (INSET_MAX - INSET_MIN)
            y = height - (top + size - (row.Display_Replicate_2_Growth_Inhibition - INSET_MIN) * size / (INSET_MAX - INSET_MIN))
            radius, handle = POINT_RADIUS, POINT_RADIUS * 0.552284749831
            commands.append(f"{x+radius:.6f} {y:.6f} m {x+radius:.6f} {y+handle:.6f} {x+handle:.6f} {y+radius:.6f} {x:.6f} {y+radius:.6f} c {x-handle:.6f} {y+radius:.6f} {x-radius:.6f} {y+handle:.6f} {x-radius:.6f} {y:.6f} c {x-radius:.6f} {y-handle:.6f} {x-handle:.6f} {y-radius:.6f} {x:.6f} {y-radius:.6f} c {x+handle:.6f} {y-radius:.6f} {x+radius:.6f} {y-handle:.6f} {x+radius:.6f} {y:.6f} c f")
        commands.extend(("Q /GS0 gs 0 0 0 1 scn 0 0 0 1 SCN 0.45 w [] 0 d",))
        rectangle(left, top, size, size, "S")
        for value in (75, 90, 105):
            x = left + (value - INSET_MIN) * size / (INSET_MAX - INSET_MIN)
            y = top + size - (value - INSET_MIN) * size / (INSET_MAX - INSET_MIN)
            line(x, top + size, x, top + size + 1.4)
            line(left - 1.4, y, left, y)
            text(x, top + size + 6, value)
            text(left - 2.8, y + 1.6, value, align="right")
    commands.append("Q")
    return ("\n".join(commands) + "\n").encode()
