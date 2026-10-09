"""Export a ROC panel using the Figure 2 vector template."""
from io import BytesIO
from pathlib import Path
import re
import pymupdf
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject
from shared.figure2_rendering import update_roc
from shared.rendering import compact_svg

METHODS = Path(__file__).resolve().parents[1]


def render_roc_panel(curves, output_stem, gi60_only=False):
    page = PdfReader(METHODS / 'shared/assets/figure2/figure_2_vector_layout.pdf').pages[0]
    content = update_roc(page.get_contents().get_data(), curves)
    if gi60_only:
        color = b'0.357 0.816 0 0  SCN'
        start = content.index(color, content.index(b'0.883 0.43 0.418 0.121  SCN'))
        end = content.index(b'S\nQ\n', start) + len(b'S\nQ\n')
        content = content[:start] + content[end:]
        legend = b'0.357 0.816 0 0  SCN\nq 1 0 0 1 63.7597 465.2715 cm\n0 0 m\n4.536 0 l\n9.072 0 l\nS\nQ\nq\n0 g\n0 G\n0 TL/Fm17 Do\nQ\n'
        if content.count(legend) != 1:
            raise ValueError('Expected one GI80 legend anchor.')
        content = content.replace(legend, b'')
    visited = set()

    def update_forms(resources):
        for name, reference in resources.get('/XObject', {}).items():
            form = reference.get_object()
            if id(form) in visited or form.get('/Subtype') != '/Form':
                continue
            visited.add(id(form))
            data = form.get_data()
            if name in ('/Fm16', '/Fm17'):
                model = 'GI60' if name == '/Fm16' else 'GI80'
                data = re.sub(rb'AUC = 0\.\d{3}', f"AUC = {curves[model]['auc']:.3f}".encode(), data)
            objects = form.get('/Resources', {}).get('/XObject', {})
            for object_name, object_reference in list(objects.items()):
                if object_reference.get_object().get('/Subtype') == '/Image':
                    data = re.sub(re.escape(object_name.encode()) + rb'\s+Do\b', b'', data)
                    del objects[object_name]
            form.set_data(data)
            if '/Resources' in form:
                update_forms(form['/Resources'])

    update_forms(page['/Resources'])
    stream = DecodedStreamObject()
    stream.set_data(content)
    page[NameObject('/Contents')] = stream.flate_encode()
    writer = PdfWriter()
    writer.add_page(page)
    buffer = BytesIO()
    writer.write(buffer)
    document = pymupdf.open(stream=buffer.getvalue(), filetype='pdf')
    page = document[0]
    page.set_cropbox(pymupdf.Rect(0, 0, 158, 173))
    output_stem = Path(output_stem).resolve()
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    document.save(output_stem.with_suffix('.pdf'), garbage=3, deflate=True)
    output_stem.with_suffix('.svg').write_text(compact_svg(page.get_svg_image(text_as_path=True)), encoding='utf-8')
    page.get_pixmap(dpi=240, alpha=False).save(output_stem.with_suffix('.png'))
    document.close()
