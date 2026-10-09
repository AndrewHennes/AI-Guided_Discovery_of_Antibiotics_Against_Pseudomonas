"""Reproduce Figure 1g v025 from current 8K activity labels, without CLI options."""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
OUT = ROOT / 'panel_g_all_hits_v025'
SOURCE = Path('/Users/asselism/Desktop/Collins_Lab/Repositories/AI-Guided_Antipseudomonal_Discovery/AI-Guided_Antipseudomonal_Discovery/datasets/8K_Dataset/Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv')
os.environ.setdefault('MPLCONFIGDIR', '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/matplotlib_cache')
sys.path.insert(0, str(ROOT / 'methods'))
from amine_descriptor import amine_counts, verify_examples
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

STRAINS = ['WT', 'Delta40', 'Delta3', 'Delta6']
COUNTS = dict(zip(STRAINS, [8, 5, 39, 124]))
LABELS = dict(zip(STRAINS, ['WT', 'Δ40', 'Δ3', 'Δ6']))
RULES = {'WT': 'WT = 1', 'Delta40': 'Delta40 = 1',
         'Delta3': 'Delta3 = 1',
         'Delta6': 'Delta6 = 1'}
COLORS = {'WT': (.07547, .65324, .31118), 'Delta40': (.97905, .63760, .11945),
          'Delta3': (.75422, .23545, .17137), 'Delta6': (.66175, .37743, .65107)}
PROPERTIES = [('TPSA', 'TPSA (Å²)', (-8, 220), [0, 50, 100, 150, 200]),
              ('Amine_count', 'Number of amines', (-.2, 3.25), [0, 1, 2, 3]),
              ('logP', 'LogP', (0, 8), [0, 2, 4, 6, 8])]
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6, 'axes.linewidth': .55,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'svg.hashsalt': 'figure1g_all_hits_v025',
                     'savefig.facecolor': 'white', 'axes.unicode_minus': False})


def canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def make_data():
    columns = {s: f'{s}_Binarized_at_80_Percent_Growth_Inhibition' for s in STRAINS}
    raw = pd.read_csv(SOURCE, usecols=['Standardized_SMILES', 'Catalog ID', *columns.values()], low_memory=False)
    raw.insert(0, 'Input_Row_Index', np.arange(len(raw)))
    assert raw[list(columns.values())].isin([0, 1]).all().all()
    raw['Canonical_SMILES'] = raw.Standardized_SMILES.map(canonical)
    grouped = raw.groupby('Canonical_SMILES', sort=True)
    calls = grouped[list(columns.values())].max().rename(columns={v: k for k, v in columns.items()})
    masks = {'WT': calls.WT.eq(1), 'Delta40': calls.Delta40.eq(1),
             'Delta3': calls.Delta3.eq(1),
             'Delta6': calls.Delta6.eq(1)}
    # Independent membership check uses sets of positive source rows.
    positives = {s: set(raw.loc[raw[col].eq(1), 'Canonical_SMILES']) for s, col in columns.items()}
    expected = {'WT': positives['WT'], 'Delta40': positives['Delta40'],
                'Delta3': positives['Delta3'],
                'Delta6': positives['Delta6']}
    assert positives['WT'] <= positives['Delta3'] <= positives['Delta6']
    for s in STRAINS:
        assert set(calls.index[masks[s]]) == expected[s]
        assert len(expected[s]) == COUNTS[s], (s, len(expected[s]))
    records = []
    for key in sorted(set.union(*expected.values())):
        mol = Chem.MolFromSmiles(key)
        rows = grouped.get_group(key)
        records.append({'Canonical_SMILES': key,
                        'Input_row_indices': ';'.join(rows.Input_Row_Index.astype(str)),
                        'Catalog_IDs': ';'.join(rows['Catalog ID'].astype(str)),
                        **{f'{s}_hit': int(calls.loc[key, s]) for s in STRAINS},
                        'TPSA': rdMolDescriptors.CalcTPSA(mol, includeSandP=False),
                        'logP': Crippen.MolLogP(mol), **amine_counts(mol)})
    wide = pd.DataFrame(records)
    parts = []
    for s in STRAINS:
        part = wide.loc[wide.Canonical_SMILES.isin(expected[s])].copy()
        part.insert(0, 'Group', s)
        part.insert(1, 'Group_rule', RULES[s])
        parts.append(part)
    long = pd.concat(parts, ignore_index=True)
    assert len(wide) == 125 and len(long) == 176
    assert not long.duplicated(['Group', 'Canonical_SMILES']).any()
    previous = pd.read_csv(ROOT / 'tables/8K_WT_and_Delta6_only_hit_properties.csv')
    joined = wide.merge(previous, on='Canonical_SMILES', suffixes=('', '_verified'))
    assert len(joined) == 124
    for prop in ['TPSA', 'logP', 'Amine_count', 'Primary_amines', 'Secondary_amines', 'Tertiary_amines']:
        assert np.allclose(joined[prop], joined[f'{prop}_verified'], rtol=1e-12, atol=1e-12), prop
    overlap = pd.DataFrame({a: {b: len(expected[a] & expected[b]) for b in STRAINS} for a in STRAINS})
    assert overlap.loc['WT', 'Delta40'] == 4
    wide.to_csv(OUT / 'all_group_properties.csv', index=False)
    long.to_csv(OUT / 'plotted_compounds.csv', index=False)
    pd.DataFrame([{'Group': s, 'N': COUNTS[s], 'Rule': RULES[s]} for s in STRAINS]).to_csv(OUT / 'group_counts.csv', index=False)
    overlap.to_csv(OUT / 'group_overlap.csv', index_label='Group')
    return raw, calls, wide, long


def offsets(group, keys):
    return np.array([(int.from_bytes(hashlib.sha256(f'{group}|{key}|property-points'.encode()).digest()[:8], 'big') / (2**64 - 1) - .5) * .42 for key in keys])


def draw(ax, data, spec):
    prop, title, limits, ticks = spec
    audit = []
    for position, group in enumerate(STRAINS):
        sub = data.loc[data.Group.eq(group)].sort_values('Canonical_SMILES')
        y = sub[prop].to_numpy()
        boxes = ax.boxplot([y], positions=[position], widths=.56, patch_artist=True,
                           whis=1.5, showfliers=False, manage_ticks=False,
                           boxprops={'linewidth': .45, 'edgecolor': '#303030', 'facecolor': COLORS[group], 'alpha': .28},
                           medianprops={'linewidth': .6, 'color': '#202020'},
                           whiskerprops={'linewidth': .45, 'color': '#303030'},
                           capprops={'linewidth': .45, 'color': '#303030'})
        for median in boxes['medians']:
            median.set_zorder(5)
        dots = ax.scatter(position + offsets(group, sub.Canonical_SMILES), y, s=2.2,
                          facecolors=[COLORS[group]], edgecolors='#202020', linewidths=.15, alpha=.7, zorder=4)
        dots.set_gid(f'{prop}_all_points_{group}')
        assert len(dots.get_offsets()) == COUNTS[group]
        assert np.isfinite(y).all() and (y > limits[0]).all() and (y < limits[1]).all()
        q1, median, q3 = np.quantile(y, [.25, .5, .75])
        within = y[(y >= q1 - 1.5 * (q3 - q1)) & (y <= q3 + 1.5 * (q3 - q1))]
        audit.append({'Property': prop, 'Group': group, 'N': len(y), 'Points_drawn': len(y),
                      'Q1': q1, 'Median': median, 'Q3': q3, 'Lower_whisker': within.min(),
                      'Upper_whisker': within.max(), 'Minimum': y.min(), 'Maximum': y.max(),
                      'Axis_minimum': limits[0], 'Axis_maximum': limits[1], 'Points_outside_axis': 0})
    ax.set(xlim=(-.55, 3.55), ylim=limits, yticks=ticks)
    ax.set_xticks(range(4), [f'{LABELS[s]}\nn={COUNTS[s]}' for s in STRAINS])
    ax.set_title(title, fontsize=6.7, pad=5)
    ax.tick_params(axis='x', labelsize=5.4, pad=3, length=0)
    ax.tick_params(axis='y', labelsize=6, pad=2, length=2, width=.5)
    ax.spines[['top', 'right']].set_visible(False)
    return audit


def main():
    OUT.mkdir(exist_ok=True)
    n_examples = verify_examples()
    raw, calls, wide, long = make_data()
    fig = plt.figure(figsize=(257 / 72, 164 / 72))
    fig.text(.01, .96, 'Physicochemical properties', fontsize=8, weight='bold', va='top')
    audit = []
    for i, spec in enumerate(PROPERTIES):
        ax = fig.add_axes([.075 + i * .323, .23, .242, .52])
        audit.extend(draw(ax, long, spec))
    for suffix in ['pdf', 'svg', 'png']:
        metadata = {'Creator': 'Figure 1g all active compounds per strain'} if suffix == 'pdf' else {'Date': None} if suffix == 'svg' else None
        fig.savefig(OUT / f'g_TPSA_amines_LogP.{suffix}', dpi=300, metadata=metadata)
    plt.close(fig)
    audit = pd.DataFrame(audit)
    assert audit.Points_drawn.sum() == 528 and audit.Points_outside_axis.sum() == 0
    audit.to_csv(OUT / 'boxplot_and_point_audit.csv', index=False)
    report = {'source': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'source_rows': len(raw), 'unique_source_structures': len(calls),
              'group_counts': COUNTS, 'group_rules': RULES, 'unique_plotted_structures': len(wide),
              'observations_per_property': len(long), 'total_points': 528,
              'group_membership_independently_checked': True,
              'descriptor_values_checked_against_verified_table': 124,
              'amine_functional_group_examples_checked': n_examples,
              'axis_limits': {s[0]: s[2] for s in PROPERTIES}, 'points_outside_axes': 0,
              'colors_RGB': COLORS, 'versions': {'rdkit': rdkit.__version__, 'matplotlib': matplotlib.__version__, 'pandas': pd.__version__}}
    (OUT / 'validation.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
