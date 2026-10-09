from pathlib import Path
from .standalone_svg import standalone_svg
from .unmasked_svg import unmasked_svg

def finalize_svg(path):
    path = Path(path).resolve()
    path.write_text(unmasked_svg(standalone_svg(path.read_text(encoding="utf-8"))), encoding="utf-8")
