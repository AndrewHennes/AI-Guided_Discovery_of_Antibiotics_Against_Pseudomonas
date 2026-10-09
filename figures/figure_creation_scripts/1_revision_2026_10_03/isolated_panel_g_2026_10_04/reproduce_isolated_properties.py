"""Create the three requested standalone property SVGs without editing Figure 1.

Absolute inputs and outputs make this script independent of the working directory.
Each strain contributes all of its binary-labelled hits. Groups may overlap.
"""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/isolated_panel_g_2026_10_04')
SOURCE = Path('/Users/asselism/Desktop/Collins_Lab/Repositories/AI-Guided_Antipseudomonal_Discovery/AI-Guided_Antipseudomonal_Discovery/datasets/8K_Dataset/Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv')
TABLES = ROOT / 'tables'
for directory in [ROOT, TABLES]:
    directory.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/matplotlib_cache')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen, GraphDescriptors, rdMolDescriptors

STRAINS = ['WT', 'Delta40', 'Delta3', 'Delta6']
LABELS = {'WT': 'WT', 'Delta40': 'Δ40', 'Delta3': 'Δ3', 'Delta6': 'Δ6'}
COLORS = {'WT': (0.07547, 0.65324, 0.31118), 'Delta40': (0.97905, 0.63760, 0.11945),
          'Delta3': (0.75422, 0.23545, 0.17137), 'Delta6': (0.66175, 0.37743, 0.65107)}
COLUMNS = {s: f'{s}_Binarized_at_80_Percent_Growth_Inhibition' for s in STRAINS}
EXPECTED_COUNTS = {'WT': 8, 'Delta40': 5, 'Delta3': 39, 'Delta6': 124}
PROPERTIES = [
    ('Log10_IPC', r'log$_{10}$ IPC', 'g_log10_IPC', (1.5, 9.0), [2, 4, 6, 8]),
    ('logP', 'logP', 'g_logP', (0, 8.0), [0, 2, 4, 6, 8]),
    ('Aromatic_rings', 'Number of aromatic rings', 'g_aromatic_rings', (-.4, 5.4), [0, 1, 2, 3, 4, 5]),
]
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': .65,
                     'xtick.major.width': .65, 'ytick.major.width': .65,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'svg.hashsalt': 'isolated_property_panels_20261004',
                     'savefig.facecolor': 'white', 'axes.unicode_minus': True})


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def build_cohorts():
    frame = pd.read_csv(SOURCE, usecols=['Standardized_SMILES', 'Catalog ID', *COLUMNS.values()], low_memory=False)
    frame.insert(0, 'Input_Row_Index', np.arange(len(frame)))
    assert frame[list(COLUMNS.values())].notna().all().all()
    assert frame[list(COLUMNS.values())].isin([0, 1]).all().all()
    frame['Canonical_SMILES'] = frame.Standardized_SMILES.map(canonicalize)
    grouped = frame.groupby('Canonical_SMILES', sort=True)
    calls = grouped[list(COLUMNS.values())].max().rename(columns={v: k for k, v in COLUMNS.items()})
    hits = calls.loc[calls.eq(1).any(axis=1)].copy()
    molecules = []
    for key, hit_row in hits.iterrows():
        mol = Chem.MolFromSmiles(key)
        ipc = GraphDescriptors.Ipc(mol, avg=False)
        assert np.isfinite(ipc) and ipc > 0, key
        rows = grouped.get_group(key)
        record = {'Canonical_SMILES': key,
                  'Input_row_indices': ';'.join(rows.Input_Row_Index.astype(str)),
                  'Catalog_IDs': ';'.join(rows['Catalog ID'].astype(str)),
                  **{f'{s}_hit': int(hit_row[s]) for s in STRAINS},
                  'IPC': ipc, 'Log10_IPC': float(np.log10(ipc)),
                  'logP': Crippen.MolLogP(mol),
                  'Aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol)}
        molecules.append(record)
    wide = pd.DataFrame(molecules)
    groups = []
    for strain in STRAINS:
        subset = wide.loc[wide[f'{strain}_hit'].eq(1)].copy()
        subset.insert(0, 'Strain', strain)
        assert len(subset) == EXPECTED_COUNTS[strain]
        assert not subset.Canonical_SMILES.duplicated().any()
        groups.append(subset)
    long = pd.concat(groups, ignore_index=True)
    wide.to_csv(TABLES / 'all_hit_properties.csv', index=False)
    long.to_csv(TABLES / 'plotted_compounds.csv', index=False)
    counts = []
    for strain in STRAINS:
        positive = calls[strain].eq(1)
        counts.append({'Strain': strain, 'Plotted_hits': int(positive.sum()),
                       'Hits_shared_with_WT': int((positive & calls.WT.eq(1)).sum()),
                       'Hits_inactive_in_WT': int((positive & calls.WT.eq(0)).sum())})
    pd.DataFrame(counts).to_csv(TABLES / 'group_counts.csv', index=False)
    overlap = pd.DataFrame({a: {b: int((calls[a].eq(1) & calls[b].eq(1)).sum()) for b in STRAINS} for a in STRAINS})
    overlap.to_csv(TABLES / 'hit_group_overlap.csv', index_label='Strain')
    assert counts[-1]['Hits_inactive_in_WT'] == 116
    return frame, calls, wide, long


def point_offsets(strain, keys):
    # Hash-based horizontal jitter is repeatable and identical across properties.
    offsets = []
    for key in keys:
        digest = hashlib.sha256(f'{strain}|{key}|property-points'.encode()).digest()
        fraction = int.from_bytes(digest[:8], 'big') / (2 ** 64 - 1)
        offsets.append((fraction - .5) * .42)
    return np.asarray(offsets)


def draw_property(ax, long, property_spec):
    prop, title, stem, limits, ticks = property_spec
    audit = []
    for position, strain in enumerate(STRAINS):
        selected = long.loc[long.Strain.eq(strain)].sort_values('Canonical_SMILES')
        values = selected[prop].to_numpy()
        color = COLORS[strain]
        # Every observation is drawn by scatter, so boxplot fliers are suppressed
        # to avoid drawing outlying compounds twice. All groups use this rule.
        boxes = ax.boxplot([values], positions=[position], widths=.56, patch_artist=True,
                           whis=1.5, showfliers=False, manage_ticks=False,
                           boxprops={'linewidth': .7, 'edgecolor': '#303030', 'facecolor': color, 'alpha': .28},
                           medianprops={'linewidth': 1.0, 'color': '#202020'},
                           whiskerprops={'linewidth': .7, 'color': '#303030'},
                           capprops={'linewidth': .7, 'color': '#303030'})
        for median in boxes['medians']:
            median.set_zorder(5)
        x = position + point_offsets(strain, selected.Canonical_SMILES)
        dots = ax.scatter(x, values, s=7, facecolors=[color], edgecolors='#202020',
                          linewidths=.22, alpha=.70, zorder=4)
        dots.set_gid(f'{stem}_all_points_{strain}')
        assert len(dots.get_offsets()) == EXPECTED_COUNTS[strain]
        assert np.all(values >= limits[0]) and np.all(values <= limits[1]), (prop, strain)
        q1, median, q3 = np.quantile(values, [.25, .5, .75])
        within = values[(values >= q1 - 1.5 * (q3 - q1)) & (values <= q3 + 1.5 * (q3 - q1))]
        audit.append({'Property': prop, 'Strain': strain, 'N': len(values), 'Points_drawn': len(values),
                      'Q1': q1, 'Median': median, 'Q3': q3, 'Lower_whisker': within.min(),
                      'Upper_whisker': within.max(), 'Minimum': values.min(), 'Maximum': values.max(),
                      'Axis_minimum': limits[0], 'Axis_maximum': limits[1], 'Points_outside_axis': 0})
    ax.set_xlim(-.55, 3.55)
    ax.set_ylim(limits)
    ax.set_yticks(ticks)
    ax.set_xticks(range(4), [f'{LABELS[s]}\nn = {EXPECTED_COUNTS[s]}' for s in STRAINS])
    ax.set_title(title, fontsize=10, weight='normal', pad=10)
    ax.tick_params(axis='x', length=0, pad=5, labelsize=8)
    ax.tick_params(axis='y', labelsize=8, pad=3)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel('Strain hit group', fontsize=8, labelpad=8)
    return audit


def main():
    frame, calls, wide, long = build_cohorts()
    audit = []
    for spec in PROPERTIES:
        fig = plt.figure(figsize=(2.65, 2.65))
        ax = fig.add_axes([.17, .26, .79, .58])
        audit.extend(draw_property(ax, long, spec))
        fig.savefig(ROOT / f'{spec[2]}.svg', metadata={'Date': None, 'Creator': 'Figure 1 isolated property panels'})
        plt.close(fig)
    pd.DataFrame(audit).to_csv(TABLES / 'boxplot_and_point_audit.csv', index=False)
    # Combined raster preview contains exactly the three property panels.
    fig = plt.figure(figsize=(7.95, 2.65))
    for i, spec in enumerate(PROPERTIES):
        ax = fig.add_axes([(i + .17) / 3, .26, .79 / 3, .58])
        draw_property(ax, long, spec)
    fig.savefig(ROOT / 'three_property_panels_preview.png', dpi=220)
    plt.close(fig)
    report = {'source_csv': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'source_records': len(frame), 'unique_source_structures': len(calls),
              'unique_hit_structures_any_strain': len(wide), 'group_counts': EXPECTED_COUNTS,
              'group_rule': 'All binary-labelled hits in each strain; canonical structures counted once per group; groups can overlap.',
              'Delta6_hits_inactive_in_WT': 116, 'inactive_group_included': False,
              'properties': [p[0] for p in PROPERTIES], 'IPC_transformation': 'log10(GraphDescriptors.Ipc(mol, avg=False))',
              'points_displayed_for_every_group': True, 'points_sampled_or_imputed': False,
              'points_outside_axes': 0, 'duplicate_boxplot_outlier_markers': False,
              'independent_group_p_values_calculated': False,
              'standalone_svgs': [f'{p[2]}.svg' for p in PROPERTIES],
              'main_illustrator_figure_edited': False, 'colors_RGB': COLORS,
              'versions': {'python': sys.version, 'rdkit': rdkit.__version__, 'numpy': np.__version__,
                           'pandas': pd.__version__, 'matplotlib': matplotlib.__version__}}
    (TABLES / 'validation.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
