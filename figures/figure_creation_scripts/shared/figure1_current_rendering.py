"""Current Figure 1 export, with preserved native panel f and vector SVGs."""
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
import pymupdf as fitz
import figure1_current_config as config
from shared.figure1_current_analysis import prepare_data
from shared.figure1_current_svg import compact_svg


def assemble():
    doc = fitz.open(config.REFERENCE_PDF)
    reference = fitz.open(config.REFERENCE_PDF)
    page = doc[0]
    replacements = [((430, 0, 655, 123), config.DATA_DIRECTORY / 'd/d_updated_compound_counts.pdf'),
                    ((618, 125, 875, 289), config.DATA_DIRECTORY / 'g/g_TPSA_amines_LogP.pdf')]
    for bounds, _ in replacements:
        page.draw_rect(fitz.Rect(bounds), color=None, fill=(1, 1, 1), overlay=True)
    # White replacement pages preserve the surrounding native artwork.
    # The standalone SVG pass removes fully covered original graphics.
    for bounds, path in replacements:
        plot = fitz.open(path)
        page.show_pdf_page(fitz.Rect(bounds), plot, 0)
    # The panel-d letter slightly overlaps the plot's white page rectangle.
    letter = fitz.Rect(423, 0, 434, 19)
    page.show_pdf_page(letter, reference, 0, clip=letter)
    preserved = {}
    for name, bounds in {'panel_f': (423, 121, 607, 288),
                         'left_panels': (0, 0, 423, 539.274),
                         'chemical_space_maps': (423, 289, 875.015, 539.274)}.items():
        a = reference[0].get_pixmap(matrix=fitz.Matrix(2, 2), clip=fitz.Rect(bounds))
        b = page.get_pixmap(matrix=fitz.Matrix(2, 2), clip=fitz.Rect(bounds))
        preserved[name] = a.samples == b.samples
    assert all(preserved.values()), preserved
    doc.set_metadata({'title': 'Figure 1, current all-hit comparison',
                      'subject': 'Current panels a–i. Panel f and original map coordinates preserved from v025.'})
    config.OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    path = config.OUTPUT_DIRECTORY / f'{config.FIGURE_NAME}.pdf'
    doc.save(path, garbage=1, deflate=True)
    doc.close()
    reference.close()
    print('Assembled current figure; retained regions are pixel-identical.', flush=True)
    return path, preserved


def validate_svg(svg):
    root = ET.fromstring(svg)
    elements = list(root.iter())
    ids = {e.get('id') for e in elements if e.get('id')}
    assert len(ids) == sum(bool(e.get('id')) for e in elements)
    for node in elements:
        for key, value in node.attrib.items():
            if key.endswith('href'):
                assert value.startswith('#') or value.startswith('data:'), value
                if value.startswith('#'):
                    assert value[1:] in ids, value
    tags = [e.tag.rsplit('}', 1)[-1] for e in elements]
    assert 'text' not in tags
    return {'elements': len(elements), 'paths': tags.count('path'), 'uses': tags.count('use'),
            'embedded_images': tags.count('image'), 'masks': tags.count('mask'), 'clipping_paths': tags.count('clipPath'),
            'external_references': 0, 'text_as_vector_outlines': True}


def export_panel(assembled, panel):
    folder = config.METHODS_DIRECTORY / f'1{panel}'
    folder.mkdir(exist_ok=True)
    stem = folder / config.PANEL_NAMES[panel]
    doc = fitz.open(assembled)
    page = doc[0]
    page.set_cropbox(fitz.Rect(config.PANEL_BOXES[panel]) & page.rect)
    # Adjacent native panels slightly overlap their rectangular crop boundaries.
    # Clear only neighboring artwork in whitespace of the standalone export.
    cleanup = {'c': [(87, 7, 92, 29), (87, 175, 92, 196)],
               'f': [(0, 19, 6, page.rect.height)],
               'h': [(0, 19, 6, page.rect.height), (227, 43, page.rect.width, page.rect.height)],
               'i': [(0, 0, 27, 43)]}
    for bounds in cleanup.get(panel, []):
        page.draw_rect(fitz.Rect(bounds), color=None, fill=(1, 1, 1), overlay=True)
    svg = compact_svg(page.get_svg_image(text_as_path=True))
    stem.with_suffix('.svg').write_text(svg)
    report = validate_svg(svg)
    # A cleaned vector SVG creates a genuinely standalone PDF without carrying
    # invisible artwork and resources from the rest of the full figure.
    if report['clipping_paths'] or report['masks']:
        # Cairo preserves SVG dash patterns on the retained density artwork.
        import cairosvg
        cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(stem.with_suffix('.pdf')))
    else:
        vector = fitz.open(stream=svg.encode(), filetype='svg')
        standalone = fitz.open(stream=vector.convert_to_pdf(), filetype='pdf')
        standalone.save(stem.with_suffix('.pdf'), garbage=3, deflate=True)
    page.get_pixmap(dpi=config.PNG_DPI, alpha=False).save(stem.with_suffix('.png'))
    config.SVG_DIRECTORY.mkdir(parents=True, exist_ok=True)
    published = config.SVG_DIRECTORY / config.SVG_NAMES[panel]
    shutil.copy2(stem.with_suffix('.svg'), published)
    report.update(panel=panel, publication_file=str(published),
                  svg_sha256=hashlib.sha256(published.read_bytes()).hexdigest())
    assert stem.with_suffix('.svg').read_bytes() == published.read_bytes()
    doc.close()
    print(f'Exported Figure 1{panel}: {config.SVG_NAMES[panel]}', flush=True)
    return report


def export_full(assembled):
    doc = fitz.open(assembled)
    stem = config.OUTPUT_DIRECTORY / config.FIGURE_NAME
    svg = compact_svg(doc[0].get_svg_image(text_as_path=True))
    stem.with_suffix('.svg').write_text(svg)
    report = validate_svg(svg)
    doc[0].get_pixmap(dpi=config.PNG_DPI, alpha=False).save(stem.with_suffix('.png'))
    dest = config.SVG_DIRECTORY / 'Figure 1 (Complete figure).svg'
    shutil.copy2(stem.with_suffix('.svg'), dest)
    report.update(publication_file=str(dest), svg_sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    doc.close()
    return report


def render_all():
    analysis = prepare_data()
    assembled, preserved = assemble()
    panels = [export_panel(assembled, p) for p in config.PANEL_NAMES]
    report = {'layout_reference': 'Figure1_v025_all_hits_TPSA_amines_LogP.ai',
              'panel_order': list(config.PANEL_NAMES), 'preserved_artwork_pixel_identical': preserved,
              'group_counts': analysis['panel_g']['group_counts'],
              'panel_f_density_table_matches_reference': analysis['panel_f']['density_table_matches_v025'],
              'panel_g_points': 528, 'standalone_panels': panels, 'full_figure': export_full(assembled)}
    (config.METHODS_DIRECTORY / 'figure_1_current_validation.json').write_text(json.dumps(report, indent=2))
    print('Figure 1 scripts and publication SVGs are synchronized.', flush=True)
    return report


def render_panel(panel):
    if panel not in config.PANEL_NAMES:
        raise ValueError(panel)
    prepare_data()
    assembled, _ = assemble()
    return export_panel(assembled, panel)
