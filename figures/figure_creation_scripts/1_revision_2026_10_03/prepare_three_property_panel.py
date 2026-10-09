"""Fit the verified three standalone property plots into the existing panel g."""
from pathlib import Path
import json
import sys
import pandas as pd

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
ISOLATED = ROOT / 'isolated_panel_g_2026_10_04'
sys.path.insert(0, str(ISOLATED / 'methods'))
import reproduce_isolated_properties as properties
import matplotlib.pyplot as plt
plt.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
                     'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold',
                     'mathtext.fallback': None})

data = pd.read_csv(ISOLATED / 'tables/plotted_compounds.csv')
fig = plt.figure(figsize=(257 / 72, 164 / 72))
fig.text(.01, .96, 'Physicochemical properties', fontsize=8, weight='bold', va='top')
fig.text(.01, .875, 'All hits per strain; groups overlap', fontsize=6)
audit = []
for i, spec in enumerate(properties.PROPERTIES):
    ax = fig.add_axes([.065 + i * .327, .23, .25, .52])
    audit.extend(properties.draw_property(ax, data, spec))
    title = 'Aromatic rings' if spec[0] == 'Aromatic_rings' else spec[1]
    ax.set_title(title, fontsize=6.7, pad=5)
    ax.set_xlabel('')
    ax.set_xticks(range(4), [f'{properties.LABELS[s]}\nn={properties.EXPECTED_COUNTS[s]}' for s in properties.STRAINS])
    ax.tick_params(axis='x', labelsize=5.4, pad=3, length=0)
    ax.tick_params(axis='y', labelsize=6, pad=2, length=2, width=.5)
    for spine in ax.spines.values():
        spine.set_linewidth(.55)
    for patch in ax.patches:
        patch.set_linewidth(.45)
    for line in ax.lines:
        line.set_linewidth(.6 if line.get_zorder() == 5 else .45)
    for points in ax.collections:
        points.set_sizes([2.2])
        points.set_linewidths([.15])
fig.text(.01, .04, 'All compounds shown; n = hits per strain', fontsize=5.7)
for suffix in ['pdf', 'svg', 'png']:
    metadata = {'Creator': 'Three-property panel integration'} if suffix == 'pdf' else {'Date': None} if suffix == 'svg' else None
    fig.savefig(ROOT / f'panels/g_three_property_panels.{suffix}', dpi=300, metadata=metadata)
plt.close(fig)
pd.DataFrame(audit).to_csv(ROOT / 'tables/integrated_three_property_point_audit.csv', index=False)
assert sum(row['Points_drawn'] for row in audit) == 528
assert all(row['Points_outside_axis'] == 0 for row in audit)
(ROOT / 'tables/three_property_layout.json').write_text(json.dumps({
    'panel_width_pt': 257, 'panel_height_pt': 164,
    'properties': [p[0] for p in properties.PROPERTIES],
    'group_counts': properties.EXPECTED_COUNTS,
    'all_points_shown_in_each_property': 176,
    'all_points_shown_total': 528,
    'points_style_same_across_groups': True,
    'data_source': str(ISOLATED / 'tables/plotted_compounds.csv'),
}, indent=2))
print('Prepared panel g with three properties and all 528 compound markers.')
