"""Export the saved version 020 Illustrator figure and verify its visible labels."""
from pathlib import Path
import json
import pymupdf as fitz

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
src = fitz.open(ROOT / 'Figure1_v020_properties.ai')
p = src[0]
assert len(src) == 1 and abs(p.rect.height - 539.274) < .01
text = p.get_text()
required = ['Growth-inhibition distribution', 'Physicochemical properties', 'Screen n = 7,994',
            '116', '7,861', 'IPC', 'Heavy atoms', 'Number of amines', 'Density floor']
assert all(s in text for s in required), [s for s in required if s not in text]
assert text.count('Screen n = 7,994') == 4
assert 'Formal charge' not in text and 'Rotatable bonds' not in text
# The source version remains available so every region outside panel g can be
# compared directly, including panel f and the user's other artwork.
previous = fitz.open(ROOT / 'Figure1_v018_violins.ai')
regions = {'above_g': (0, 0, 875.015, 125),
           'left_of_g': (0, 125, 618, 289),
           'below_g': (0, 289, 875.015, 539.274)}
unchanged = {}
for name, rect in regions.items():
    clip = fitz.Rect(rect)
    before = previous[0].get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip)
    after = p.get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip)
    unchanged[name] = before.samples == after.samples
assert all(unchanged.values()), unchanged
out = fitz.open()
out.insert_pdf(src, links=False, annots=False)
out.set_metadata({'title': 'Figure 1, version 020', 'subject': 'Growth inhibition violins and physicochemical property box-and-whisker plots.'})
out.save(ROOT / 'Figure1_v020_properties.pdf', garbage=1, deflate=True, use_objstms=1)
p.get_pixmap(matrix=fitz.Matrix(2, 2)).save(ROOT / 'Figure1_v020_properties.png')
pix = p.get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
review = fitz.open()
page = review.new_page(width=p.rect.width, height=p.rect.height)
page.insert_image(page.rect, stream=pix.tobytes('png'))
review.set_metadata({'title': 'Figure 1, version 020, 300 dpi review copy', 'subject': 'Raster preview. Editable artwork is in the AI and vector PDF.'})
review.save(ROOT / 'Figure1_v020_review.pdf', deflate=True)
report = {'page_count': len(src), 'width_pt': p.rect.width, 'height_pt': p.rect.height,
          'screen_annotation_count': text.count('Screen n = 7,994'), 'IPC_label_found': True,
          'formal_charge_label_removed': True, 'amine_count_and_heavy_atoms_labels_found': True, 'review_pdf_dpi': 300,
          'regions_outside_panel_g_pixel_identical': unchanged,
          'exports_derived_from_saved_native_AI': True,
          'files': {n: (ROOT / n).stat().st_size for n in ['Figure1_v020_properties.ai',
              'Figure1_v020_properties.pdf', 'Figure1_v020_review.pdf', 'Figure1_v020_properties.png']}}
(ROOT / 'tables/property_descriptor_artifact_validation.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
