"""Verify data membership and the actual point elements in each standalone SVG."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/isolated_panel_g_2026_10_04')
SOURCE = Path('/Users/asselism/Desktop/Collins_Lab/Repositories/AI-Guided_Antipseudomonal_Discovery/AI-Guided_Antipseudomonal_Discovery/datasets/8K_Dataset/Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv')
REFERENCE = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/tables/all_8k_molecular_properties.csv')
expected = {'WT': 8, 'Delta40': 5, 'Delta3': 39, 'Delta6': 124}
calls = {s: f'{s}_Binarized_at_80_Percent_Growth_Inhibition' for s in expected}
source = pd.read_csv(SOURCE, usecols=['Standardized_SMILES', *calls.values()])
long = pd.read_csv(ROOT / 'tables/plotted_compounds.csv')
wide = pd.read_csv(ROOT / 'tables/all_hit_properties.csv')
assert len(wide) == 125 and not wide.Canonical_SMILES.duplicated().any()
for strain, n in expected.items():
    positive = source.loc[source[calls[strain]].eq(1), 'Standardized_SMILES']
    keys = {Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=True) for s in positive}
    group = long.loc[long.Strain.eq(strain)]
    assert len(keys) == len(group) == n
    assert set(group.Canonical_SMILES) == keys
reference = pd.read_csv(REFERENCE).set_index('Canonical_SMILES')
for prop in ['IPC', 'Log10_IPC', 'logP', 'Aromatic_rings']:
    np.testing.assert_allclose(wide[prop], reference.loc[wide.Canonical_SMILES, prop], rtol=1e-12, atol=1e-12)
ns = {'s': 'http://www.w3.org/2000/svg'}
svg_checks = {}
files = ['g_log10_IPC.svg', 'g_logP.svg', 'g_aromatic_rings.svg']
assert sorted(p.name for p in ROOT.glob('*.svg')) == sorted(files)
for name in files:
    path = ROOT / name
    svg = ET.parse(path).getroot()
    assert not svg.findall('.//s:image', ns)
    assert svg.findall('.//s:text', ns)
    counts = {}
    for strain, n in expected.items():
        group = svg.find(f".//s:g[@id='{path.stem}_all_points_{strain}']", ns)
        assert group is not None
        counts[strain] = len(group.findall('.//s:use', ns))
        assert counts[strain] == n
    svg_checks[name] = {'vector_points_per_group': counts, 'text_editable': True, 'embedded_raster_images': 0}
audit = pd.read_csv(ROOT / 'tables/boxplot_and_point_audit.csv')
assert len(audit) == 12 and audit.Points_drawn.eq(audit.N).all()
assert audit.Points_outside_axis.eq(0).all()
report = {'group_membership_matches_positive_source_records': True,
          'descriptors_match_previously_verified_values': True,
          'all_176_group_observations_in_each_svg': True,
          'exactly_three_standalone_svg_panels': True,
          'no_inactive_group': set(long.Strain) == set(expected),
          'svg_checks': svg_checks}
(ROOT / 'tables/artifact_checks.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
