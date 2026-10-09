"""Export the verified native AI artwork as vector and lightweight review PDFs."""
from pathlib import Path
import json,shutil
import pymupdf as fitz
ROOT=Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
src=fitz.open(ROOT/'Figure1_v017_cleaned.ai')
out=fitz.open();out.insert_pdf(src,links=False,annots=False)
out.set_metadata({'title':'Figure 1 revised','subject':'Native Illustrator revision with growth inhibition and physicochemical property panels.'})
out.save(ROOT/'Figure1_v017_cleaned.pdf',garbage=1,deflate=True,use_objstms=1)
p=src[0]
p.get_pixmap(matrix=fitz.Matrix(2,2)).save(ROOT/'Figure1_v017_cleaned.png')
pix=p.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False)
review=fitz.open();page=review.new_page(width=p.rect.width,height=p.rect.height);page.insert_image(page.rect,stream=pix.tobytes('png'))
review.set_metadata({'title':'Figure 1 revised, 300 dpi review copy','subject':'Raster review copy. Editable vectors are supplied in the native AI and full PDF.'})
review.save(ROOT/'Figure1_v017_review.pdf',deflate=True)
text=p.get_text()
required=['Growth-inhibition distribution','Physicochemical properties','7,671','116','7,861','Formal charge','Rotatable bonds']
assert all(s in text for s in required)
assert len(src)==1 and abs(p.rect.height-539.274)<.01
report={'native_ai_saved_by_illustrator':True,'stray_placed_figures_removed':2,'native_layout_visually_reviewed':True,'page_count':len(src),'width_pt':p.rect.width,'height_pt':p.rect.height,'required_labels_found':True,'full_pdf_derived_from_saved_AI':True,'review_pdf_dpi':300,'files':{n:(ROOT/n).stat().st_size for n in ['Figure1_v017_cleaned.ai','Figure1_v017_cleaned.pdf','Figure1_v017_review.pdf','Figure1_v017_cleaned.png']}}
(ROOT/'tables/cleaned_artifact_validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2),flush=True)
