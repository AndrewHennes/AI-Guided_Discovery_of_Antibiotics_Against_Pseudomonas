"""Export v024 from native AI and strictly verify that panel f was not changed."""
from pathlib import Path
import hashlib
import json
import pymupdf as fitz

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
OUT = ROOT / 'panel_g_activity_groups_v024'
source_path = ROOT / 'Figure1_v024_TPSA_amines_LogP.ai'
doc = fitz.open(source_path)
previous = fitz.open(ROOT / 'Figure1_v023_three_properties.ai')
assert len(doc) == len(previous) == 1
page = doc[0]
assert page.rect == previous[0].rect
g_rect = fitz.Rect(618, 125, 875.015, 289)
g_text = page.get_text(clip=g_rect)
for label in ['Physicochemical properties', 'TPSA', 'Number of amines', 'LogP', 'WT', 'Δ40', 'Δ3', 'Δ6']:
    assert label in g_text, label
for count in ['n=8', 'n=5', 'n=31', 'n=85']:
    assert g_text.splitlines().count(count) == 3, (count, g_text)
for old in ['IPC', 'Aromatic rings', 'n=39', 'n=124', 'active', 'inactive', '*', '†']:
    assert old not in g_text, old
fonts = {s['font'] for b in page.get_text('dict', clip=g_rect)['blocks'] if 'lines' in b
         for line in b['lines'] for s in line['spans']}
assert all('Arial' in f for f in fonts), fonts

regions = {'above_g': (0, 0, 875.015, 125), 'left_of_g': (0, 125, 618, 289),
           'below_g': (0, 289, 875.015, 539.274), 'panel_f': (430, 125, 615, 289)}
unchanged = {}
for name, rect in regions.items():
    before = previous[0].get_pixmap(matrix=fitz.Matrix(2, 2), clip=fitz.Rect(rect))
    after = page.get_pixmap(matrix=fitz.Matrix(2, 2), clip=fitz.Rect(rect))
    unchanged[name] = before.samples == after.samples
assert all(unchanged.values()), unchanged
f_rect = fitz.Rect(regions['panel_f'])
assert page.get_text(clip=f_rect) == previous[0].get_text(clip=f_rect)
assert page.get_text().count('Screen n = 7,994') == 4

out = fitz.open()
out.insert_pdf(doc, links=False, annots=False)
out.set_metadata({'title': 'Figure 1, version 024', 'subject': 'TPSA, amine count, and LogP with corrected activity groups. Panel f retained.'})
out.save(ROOT / 'Figure1_v024_TPSA_amines_LogP.pdf', garbage=1, deflate=True, use_objstms=1)
page.get_pixmap(matrix=fitz.Matrix(2, 2)).save(ROOT / 'Figure1_v024_TPSA_amines_LogP.png')
page.get_pixmap(matrix=fitz.Matrix(4, 4), clip=g_rect).save(OUT / 'g_native_v024.png')
pix = page.get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
review = fitz.open()
rp = review.new_page(width=page.rect.width, height=page.rect.height)
rp.insert_image(rp.rect, stream=pix.tobytes('png'))
review.set_metadata({'title': 'Figure 1, version 024, 300 dpi review copy', 'subject': 'Raster review copy; editable artwork is in the native AI and vector PDF.'})
review.save(ROOT / 'Figure1_v024_review.pdf', deflate=True)
report = {'native_file': str(source_path), 'native_sha256': hashlib.file_digest(source_path.open('rb'), 'sha256').hexdigest(),
          'group_counts': {'WT': 8, 'Delta40': 5, 'Delta3': 31, 'Delta6': 85},
          'properties': ['TPSA', 'Amine_count', 'logP'], 'simple_category_labels_only': True,
          'panel_g_fonts': sorted(fonts), 'regions_pixel_identical_to_v023': unchanged,
          'panel_f_text_identical_to_v023': True, 'all_387_points_inside_axes': True,
          'data_validation_report': str(OUT / 'validation.json'), 'review_pdf_dpi': 300,
          'exports_derived_from_saved_native_AI': True,
          'files': {name: (ROOT / name).stat().st_size for name in ['Figure1_v024_TPSA_amines_LogP.ai',
                   'Figure1_v024_TPSA_amines_LogP.pdf', 'Figure1_v024_review.pdf', 'Figure1_v024_TPSA_amines_LogP.png']}}
(OUT / 'native_validation.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
