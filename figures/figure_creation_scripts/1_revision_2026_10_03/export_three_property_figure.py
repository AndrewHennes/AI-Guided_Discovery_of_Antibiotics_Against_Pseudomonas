"""Export and check the saved native version 023 Illustrator figure."""
from pathlib import Path
import hashlib
import json
import pymupdf as fitz

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
source_path = ROOT / 'Figure1_v023_three_properties.ai'
src = fitz.open(source_path)
page = src[0]
assert len(src) == 1 and abs(page.rect.height - 539.274) < .01
panel = fitz.Rect(618, 125, 875.015, 289)
text = page.get_text()
panel_text = page.get_text(clip=panel)
required = ['Physicochemical properties', 'IPC', 'logP', 'Aromatic rings',
            'All hits per strain; groups overlap', 'All compounds shown',
            'WT', 'Δ40', 'Δ3', 'Δ6']
assert all(label in panel_text for label in required), [label for label in required if label not in panel_text]
for count in ['n=8', 'n=5', 'n=39', 'n=124']:
    assert panel_text.count(count) == 3, (count, panel_text.count(count))
for old_label in ['Heavy atoms', 'Number of amines', 'TPSA', 'Inactive', 'Δ6-only']:
    assert old_label not in panel_text, old_label
assert text.count('Screen n = 7,994') == 4

previous = fitz.open(ROOT / 'Figure1_v022_heavy_atoms.ai')
regions = {'above_g': (0, 0, 875.015, 125),
           'left_of_g': (0, 125, 618, 289),
           'below_g': (0, 289, 875.015, 539.274)}
unchanged = {}
render_comparison = {}
for name, rect in regions.items():
    clip = fitz.Rect(rect)
    before = previous[0].get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip)
    after = page.get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip)
    unchanged[name] = before.samples == after.samples
    a, b = before.samples, after.samples
    different = [i // 3 for i in range(0, len(a), 3) if a[i:i + 3] != b[i:i + 3]]
    maximum_delta = max(abs(x - y) for x, y in zip(a, b)) if different else 0
    bounds = None
    if different:
        xs = [i % before.width for i in different]
        ys = [i // before.width for i in different]
        bounds = [min(xs) / 2 + clip.x0, min(ys) / 2 + clip.y0,
                  max(xs) / 2 + clip.x0, max(ys) / 2 + clip.y0]
    render_comparison[name] = {'different_pixels': len(different),
                               'maximum_channel_delta_out_of_255': maximum_delta,
                               'difference_bounds_pt': bounds}
assert unchanged['above_g'] and unchanged['below_g'], unchanged
# Reopening the CMYK Illustrator file and saving it in 30.8.2 causes a small
# colour conversion in the RGB violin paths. Check its measured extent and
# vector geometry explicitly instead of describing it as pixel-identical.
f_rect = fitz.Rect(430, 125, 615, 289)
left = render_comparison['left_of_g']
assert left['maximum_channel_delta_out_of_255'] <= 4, left
if left['difference_bounds_pt']:
    assert f_rect.contains(fitz.Rect(left['difference_bounds_pt'])), left
drawings = [[d for d in p.get_cdrawings()
             if d.get('rect') and f_rect.contains(fitz.Rect(d['rect']))]
            for p in [previous[0], page]]
geometry_keys = ['items', 'rect', 'type', 'width', 'lineCap', 'lineJoin',
                 'dashes', 'closePath', 'even_odd', 'fill_opacity', 'stroke_opacity']
assert len(drawings[0]) == len(drawings[1]) == 32
assert all(all(a.get(k) == b.get(k) for k in geometry_keys)
           for a, b in zip(*drawings)), 'Panel f geometry or opacity changed.'
panel_fonts = {span['font'] for block in page.get_text('dict', clip=panel)['blocks']
               if 'lines' in block for line in block['lines'] for span in line['spans']}
assert not any('DejaVu' in font for font in panel_fonts), panel_fonts

out = fitz.open()
out.insert_pdf(src, links=False, annots=False)
out.set_metadata({'title': 'Figure 1, version 023', 'subject': 'Three physicochemical properties across WT, Delta40, Delta3 and Delta6 hits.'})
out.save(ROOT / 'Figure1_v023_three_properties.pdf', garbage=1, deflate=True, use_objstms=1)
page.get_pixmap(matrix=fitz.Matrix(2, 2)).save(ROOT / 'Figure1_v023_three_properties.png')
page.get_pixmap(matrix=fitz.Matrix(4, 4), clip=panel).save(ROOT / 'panels/g_native_v023.png')
pix = page.get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
review = fitz.open()
preview_page = review.new_page(width=page.rect.width, height=page.rect.height)
preview_page.insert_image(preview_page.rect, stream=pix.tobytes('png'))
review.set_metadata({'title': 'Figure 1, version 023, 300 dpi review copy', 'subject': 'Raster review copy. Editable artwork is in the AI and vector PDF.'})
review.save(ROOT / 'Figure1_v023_review.pdf', deflate=True)
report = {
    'page_count': len(src), 'width_pt': page.rect.width, 'height_pt': page.rect.height,
    'screen_annotation_count': text.count('Screen n = 7,994'),
    'properties': ['log10 IPC', 'logP', 'Aromatic rings'],
    'group_counts': {'WT': 8, 'Delta40': 5, 'Delta3': 39, 'Delta6': 124},
    'inactive_group_removed': True,
    'point_audit': 'tables/integrated_three_property_point_audit.csv',
    'regions_outside_panel_g_pixel_identical': unchanged,
    'render_comparison_outside_panel_g': render_comparison,
    'panel_f_vector_geometry_and_opacity_identical': True,
    'panel_g_fonts': sorted(panel_fonts),
    'exports_derived_from_saved_native_AI': True, 'review_pdf_dpi': 300,
    'native_ai_sha256': hashlib.file_digest(source_path.open('rb'), 'sha256').hexdigest(),
    'files': {name: (ROOT / name).stat().st_size for name in [
        'Figure1_v023_three_properties.ai', 'Figure1_v023_three_properties.pdf',
        'Figure1_v023_review.pdf', 'Figure1_v023_three_properties.png']},
}
(ROOT / 'tables/three_property_native_validation.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
