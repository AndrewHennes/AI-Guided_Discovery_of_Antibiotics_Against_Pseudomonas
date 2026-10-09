"""Vector PDF preview assembled from the supplied AI file's PDF representation."""
from pathlib import Path
import pymupdf as fitz
ROOT=Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
SOURCE=Path('/Users/asselism/Downloads/Figure1_v014 (3).ai')
SHIFT=70
src=fitz.open(SOURCE);out=fitz.open();p=out.new_page(width=src[0].rect.width,height=src[0].rect.height+SHIFT)
# Preserve all original artwork except the count panel and map placement.
for r in [fitz.Rect(0,0,433,469.274),fitz.Rect(650,0,875.015,221)]:
    p.show_pdf_page(r,src,0,clip=r)
# Preserve the original two chemical-space maps and their legends as vectors.
r=fitz.Rect(433,221,875.015,469.274)
p.show_pdf_page(fitz.Rect(r.x0,r.y0+SHIFT,r.x1,r.y1+SHIFT),src,0,clip=r)
for stem,r in [('d_updated_compound_counts',fitz.Rect(430,0,650,123)),('f_inhibition_ridgelines',fitz.Rect(430,125,615,289)),('g_physicochemical_properties',fitz.Rect(618,125,875,289))]:
    with fitz.open(ROOT/'panels'/f'{stem}.pdf') as panel:p.show_pdf_page(r,panel,0)
# Panel labels use the source label font appearance.
for t,x,y in [('f',424,138),('g',607,138),('h',424,305),('i',677,305)]:
    p.insert_text((x,y),t,fontname='hebo',fontsize=13.5,color=(.106,.110,.094))
out.save(ROOT/'Figure1_v016_revised.pdf',garbage=1,deflate=True,use_objstms=1)
p.get_pixmap(matrix=fitz.Matrix(2,2)).save(ROOT/'Figure1_v016_revised.png')
print(ROOT/'Figure1_v016_revised.pdf')
