"""Preview panel placement from the saved artwork before the native Illustrator save."""
from pathlib import Path
import pymupdf as fitz

ROOT=Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
source=fitz.open(ROOT/'Figure1_v022_heavy_atoms.ai')
panel=fitz.open(ROOT/'panels/g_three_property_panels.pdf')
page=source[0]
target=fitz.Rect(618,125,875,289)
page.draw_rect(target,color=None,fill=(1,1,1),overlay=True)
page.show_pdf_page(target,panel,0,overlay=True)
page.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False).save(ROOT/'Figure1_v023_layout_preview.png')
preview=fitz.open()
out=preview.new_page(width=page.rect.width,height=page.rect.height)
pix=page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False)
out.insert_image(out.rect,stream=pix.tobytes('png'))
preview.set_metadata({'title':'Figure 1, three-property layout preview','subject':'Placement preview pending the native Illustrator save.'})
preview.save(ROOT/'Figure1_v023_layout_preview.pdf',deflate=True)
print('Prepared the complete figure layout preview.')
