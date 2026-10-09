"""Validate scientific calculations and display choices for the violin revision."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import GraphDescriptors
from violin_panels import DENSITY_FLOOR, PROPERTY_LIMITS, SCREEN_N

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
t = ROOT / 'tables'
p = pd.read_csv(t / 'all_8k_molecular_properties.csv')
assert p.Category.value_counts().to_dict() == {'Inactive_in_WT_and_Delta6': 7861, 'Delta6_only_hit': 116, 'WT_hit': 8}
assert len(p) == 7985 and not p.Canonical_SMILES.duplicated().any()
assert np.isfinite(p[['IPC', 'Log10_IPC']]).all().all() and p.IPC.gt(0).all()
np.testing.assert_allclose(p.Log10_IPC, np.log10(p.IPC), rtol=1e-13)
# Verify the exported descriptor against structures from every group, including all WT hits.
sample = pd.concat([p.loc[p.Category.eq('WT_hit')], p.loc[~p.Category.eq('WT_hit')].sample(100, random_state=1826)])
for row in sample.itertuples():
    value = GraphDescriptors.Ipc(Chem.MolFromSmiles(row.Canonical_SMILES), avg=False)
    np.testing.assert_allclose(row.IPC, value, rtol=1e-12)
audit = pd.read_csv(t / 'property_violin_display_audit.csv')
for row in audit.itertuples():
    a = p.loc[p.Category.eq(row.Group), row.Property]
    lo, hi = PROPERTY_LIMITS[row.Property]
    assert row.N == len(a)
    assert row.Below_display == a.lt(lo).sum()
    assert row.Above_display == a.gt(hi).sum()
assert audit.loc[audit.Group.eq('WT_hit'), ['Below_display', 'Above_display']].sum().sum() == 0
d = pd.read_csv(t / 'inhibition_violin_densities.csv')
settings = json.loads((t / 'inhibition_violin_settings.json').read_text())
assert settings['screen_n_annotation'] == SCREEN_N == 7994
assert settings['measured_unique_structures'] == {'WT': 7671, 'Delta40': 7985, 'Delta3': 7985, 'Delta6': 7985}
integrals = {}
for strain in ['WT', 'Delta40', 'Delta3', 'Delta6']:
    density = d[f'{strain}_density'].to_numpy()
    width = d[f'{strain}_half_width'].to_numpy()
    integral = np.trapezoid(density, d.Mean_growth_inhibition_percent)
    assert abs(integral - 1) < 1e-5, (strain, integral)
    expected = .405 * np.maximum(np.log10(np.maximum(density, DENSITY_FLOOR) / DENSITY_FLOOR), 0) / np.log10(settings['shared_peak_density'] / DENSITY_FLOOR)
    np.testing.assert_allclose(width, expected, atol=1e-12)
    integrals[strain] = float(integral)
report = {'descriptor_values_checked_from_structures': len(sample), 'all_log10_IPC_values_finite': True,
          'direct_log10_no_pseudocount': True, 'WT_hits': 8, 'Delta6_only_hits': 116,
          'inactive_in_both': 7861, 'all_WT_points_inside_display': True,
          'outside_range_counts_verified': True, 'density_integrals': integrals,
          'shared_log_density_widths_verified': True, 'screen_annotation_n': SCREEN_N,
          'measured_density_counts_documented_separately': True}
(t / 'violin_analysis_checks.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
