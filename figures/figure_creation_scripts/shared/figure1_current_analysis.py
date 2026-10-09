"""Current hit counts, unchanged inhibition-density calculations, and properties."""
import hashlib
import json
import os
from functools import lru_cache
from pathlib import Path
import figure1_current_config as config
os.environ.setdefault('MPLCONFIGDIR', str(config.DATA_DIRECTORY / 'matplotlib_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from rdkit import Chem
from shared import figure1_properties as properties
from shared.figure1_inhibition import plot_inhibition_violins

STRAINS = ['WT', 'Delta40', 'Delta3', 'Delta6']
COLORS = {**properties.COLORS, 'Delta5': (.36083, .75999, .87846)}
LABELS = {'WT': 'WT', 'Delta40': 'ΔOMP40', 'Delta3': 'ΔE-pump3', 'Delta6': 'ΔE-pump6', 'Delta5': 'ΔE-pump5'}
COLUMNS = {'8K': {s: f'{s}_Binarized_at_80_Percent_Growth_Inhibition' for s in STRAINS},
           '39K': {'WT': 'PA_Hit_Binarized_at_80_Percent_Growth_Inhibition'},
           'ECBD': {'WT': 'Pseudomonas_aeruginosa_Binarized_Hit'},
           'COADD': {'WT': 'PA_Binarized_at_80_Percent_Growth_Inhibition',
                     'Delta5': 'PA_MexB_Binarized_at_80_Percent_Growth_Inhibition'}}


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
    directory = config.DATA_DIRECTORY / 'd'
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for dataset, columns in COLUMNS.items():
        frame = pd.read_csv(config.INPUT_DIRECTORY / f'{dataset}.csv.gz', low_memory=False)
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


def inhibition():
    directory = config.DATA_DIRECTORY / 'f'
    directory.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(config.INPUT_DIRECTORY / '8K.csv.gz', low_memory=False)
    frame['Canonical_SMILES'] = frame.Standardized_SMILES.map(canonical)
    columns = [f'Average_{s}_Percent_Growth_Inhibition' for s in STRAINS]
    grouped = frame.groupby('Canonical_SMILES', sort=True)
    values = grouped[columns].mean()
    expected = pd.read_csv(config.ASSET_DIRECTORY / 'reference_tables/inhibition_values_by_structure.csv', index_col=0)
    pd.testing.assert_frame_equal(values, expected, check_exact=False, rtol=1e-10, atol=1e-10)
    summary = []
    for strain, column in zip(STRAINS, columns):
        summary.append({'Strain': strain, 'Measured_unique_structures': int(values[column].notna().sum()),
                        'Binary_hit_structures': int(grouped[COLUMNS['8K'][strain]].max().eq(1).sum())})
    summary = pd.DataFrame(summary)
    values.to_csv(directory / 'inhibition_values_by_structure.csv')
    summary.to_csv(directory / 'inhibition_summary.csv', index=False)
    plot_inhibition_violins(values, summary, STRAINS, LABELS, COLORS,
        lambda fig, stem: save_plot(fig, directory, stem), directory)
    actual = pd.read_csv(directory / 'inhibition_violin_densities.csv')
    reference = pd.read_csv(config.ASSET_DIRECTORY / 'reference_tables/inhibition_violin_densities.csv')
    pd.testing.assert_frame_equal(actual, reference, check_exact=False, rtol=1e-10, atol=1e-12)
    # The published panel f is exported from the saved native artwork. Its
    # independently recalculated values and density curves must remain identical.
    return {'measured_counts': dict(zip(summary.Strain, summary.Measured_unique_structures)),
            'screen_n_annotation': 7994, 'density_table_matches_v025': True,
            'published_artwork_preserved': True}


def prepare_data():
    config.DATA_DIRECTORY.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': .6,
                         'xtick.major.width': .6, 'ytick.major.width': .6,
                         'xtick.major.size': 2, 'ytick.major.size': 2,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    count_report = counts()
    inhibition_report = inhibition()
    plt.rcParams['axes.linewidth'] = .55
    properties.main()
    report = {'dataset_hit_counts': count_report, 'panel_f': inhibition_report,
              'panel_g': json.loads((config.DATA_DIRECTORY / 'g/validation.json').read_text()),
              'input_checksums': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in config.INPUT_DIRECTORY.glob('*.csv.gz')}}
    (config.DATA_DIRECTORY / 'analysis_validation.json').write_text(json.dumps(report, indent=2))
    return report
