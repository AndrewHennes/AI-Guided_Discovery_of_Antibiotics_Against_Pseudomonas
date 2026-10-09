"""Restore the verified version 020 artwork for the requested return to heavy atoms.

This is an exact restoration of that saved figure, not a way to apply fresh
analysis results. For regenerated panels, use the normal Illustrator update
script instead. Version 021 was confirmed saved and unmodified before this
restoration; its only figure change from 020 was the alcohol-count subplot.
"""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
source = ROOT / 'Figure1_v020_properties.ai'
target = ROOT / 'Figure1_v022_heavy_atoms.ai'
assert not target.exists(), 'Refusing to overwrite an existing restored document.'
shutil.copy2(source, target)
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
assert source_hash == hashlib.sha256(target.read_bytes()).hexdigest()
(ROOT / 'tables/heavy_atom_restore_provenance.json').write_text(json.dumps({
    'source_native_artwork': str(source), 'restored_native_artwork': str(target),
    'sha256': source_hash, 'restoration_byte_identical': True,
    'current_version_021_confirmed_saved_without_intervening_edits': True,
    'current_property_values_match_verified_version_020': True,
}, indent=2))
print(target)
