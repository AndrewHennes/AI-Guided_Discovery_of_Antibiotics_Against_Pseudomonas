"""Figure 1d generating code: current unique hit counts and grouped bars.

Uses bundled source snapshots and writes to this panel's generated_from_data
folder. No CLI arguments; all paths resolve from this file's absolute location.
"""
from pathlib import Path
from types import SimpleNamespace
from functools import lru_cache
import hashlib
import json
import os

PANEL_DIRECTORY = Path(__file__).resolve().parent
METHODS_DIRECTORY = PANEL_DIRECTORY.parent
OUTPUT_DIRECTORY = PANEL_DIRECTORY / 'generated_from_data'
ASSET_DIRECTORY = METHODS_DIRECTORY / 'shared' / 'figure1_current'
INPUT_DIRECTORY = ASSET_DIRECTORY / 'inputs'
config = SimpleNamespace(ASSET_DIRECTORY=ASSET_DIRECTORY, INPUT_DIRECTORY=INPUT_DIRECTORY)
os.environ.setdefault('MPLCONFIGDIR', str(OUTPUT_DIRECTORY / 'matplotlib_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from rdkit import Chem

STRAINS = ['WT', 'Delta40', 'Delta3', 'Delta6']
COLORS = {'WT': (.07547, .65324, .31118), 'Delta40': (.97905, .63760, .11945),
          'Delta3': (.75422, .23545, .17137), 'Delta6': (.66175, .37743, .65107),
          'Delta5': (.36083, .75999, .87846)}
LABELS = {'WT': 'WT', 'Delta40': 'ΔOMP40', 'Delta3': 'ΔE-pump3', 'Delta6': 'ΔE-pump6', 'Delta5': 'ΔE-pump5'}
COLUMNS = {'8K': {s: f'{s}_Binarized_at_80_Percent_Growth_Inhibition' for s in STRAINS},
           '39K': {'WT': 'PA_Hit_Binarized_at_80_Percent_Growth_Inhibition'},
           'ECBD': {'WT': 'Pseudomonas_aeruginosa_Binarized_Hit'},
           'COADD': {'WT': 'PA_Binarized_at_80_Percent_Growth_Inhibition',
                     'Delta5': 'PA_MexB_Binarized_at_80_Percent_Growth_Inhibition'}}
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': .6,
                     'xtick.major.width': .6, 'ytick.major.width': .6,
                     'xtick.major.size': 2, 'ytick.major.size': 2,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.unicode_minus': False, 'savefig.facecolor': 'white'})

@lru_cache(None)
def canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def save_plot(fig, directory, stem):
    directory.mkdir(parents=True, exist_ok=True)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(directory / f'{stem}.{ext}', dpi=240)
    plt.close(fig)

def counts():
    directory = OUTPUT_DIRECTORY
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for dataset, columns in COLUMNS.items():
        frame = pd.read_csv(config.INPUT_DIRECTORY / f'{dataset.lower()}.csv.gz', low_memory=False)
        smiles = 'Standardized_SMILES' if 'Standardized_SMILES' in frame else 'SMILES'
        for strain, column in columns.items():
            assert frame[column].dropna().isin([0, 1]).all()
            # Unique positive structures implement OR across duplicate records.
            positive = frame.loc[frame[column].eq(1), smiles].map(canonical)
            records.append({'Dataset': dataset, 'Strain': strain,
                            'Unique_hit_structures': positive.nunique(), 'Hit_records': len(positive)})
    result = pd.DataFrame(records)
    expected = pd.read_csv(config.ASSET_DIRECTORY / 'reference_tables/updated_hit_counts.csv')
    pd.testing.assert_frame_equal(result, expected[result.columns], check_dtype=False)
    result.to_csv(directory / 'updated_hit_counts.csv', index=False)
    fig, ax = plt.subplots(figsize=(225 / 72, 123 / 72))
    fig.subplots_adjust(left=.17, right=.99, bottom=.30, top=.79)
    for i, dataset in enumerate(COLUMNS):
        subset = result.loc[result.Dataset.eq(dataset)]
        for j, row in enumerate(subset.itertuples()):
            x = i + (j - (len(subset) - 1) / 2) * .17
            ax.bar(x, row.Unique_hit_structures, .145, color=COLORS[row.Strain], edgecolor='black', linewidth=.4)
            ax.text(x, row.Unique_hit_structures + 2, str(row.Unique_hit_structures), ha='center', va='bottom', fontsize=7)
    ax.set(ylim=(0, 140), yticks=[0, 40, 80, 120])
    ax.set_ylabel('Number of unique hits', fontsize=7, labelpad=2)
    ax.set_xticks(range(4), ['Internal 8K', 'Internal 39K', 'ECBD', 'CO-ADD'], rotation=28, ha='right', fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=6, pad=1.5)
    fig.legend([Patch(facecolor=COLORS[s], edgecolor='black', linewidth=.4) for s in [*STRAINS, 'Delta5']],
               [LABELS[s] for s in [*STRAINS, 'Delta5']], loc='upper center', bbox_to_anchor=(.57, 1.01),
               ncol=3, frameon=False, fontsize=6, handlelength=1, columnspacing=.75, handletextpad=.3, labelspacing=.3)
    save_plot(fig, directory, 'd_updated_compound_counts')
    return records

if __name__ == '__main__':
    result = counts()
    (OUTPUT_DIRECTORY / 'validation.json').write_text(json.dumps({'counts_match_final_figure': True, 'counts': result}, indent=2))
    print('Figure 1d generated; all unique-hit counts match the final figure.')
