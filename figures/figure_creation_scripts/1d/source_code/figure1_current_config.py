"""Current Figure 1 (v025); all paths resolve absolutely, independent of cwd."""
from pathlib import Path

METHODS_DIRECTORY = Path(__file__).resolve().parent
ASSET_DIRECTORY = METHODS_DIRECTORY / 'shared' / 'figure1_current'
INPUT_DIRECTORY = ASSET_DIRECTORY / 'inputs'
DATA_DIRECTORY = METHODS_DIRECTORY / 'current_figure_1_data'
OUTPUT_DIRECTORY = METHODS_DIRECTORY / 'full_figure'
SVG_DIRECTORY = METHODS_DIRECTORY.parent / 'figure_svg_files'
REFERENCE_PDF = ASSET_DIRECTORY / 'figure1_v025_reference.pdf'
FIGURE_NAME = 'antipseudomonal_screening_and_chemical_space_figure_1'
PNG_DPI = 180
PANEL_NAMES = {
    'a': 'antipseudomonal_compound_screening_workflow_a',
    'b': 'replicate_growth_inhibition_correlations_across_four_strains_b',
    'c': 'replicate_growth_inhibition_zoomed_views_c',
    'd': 'unique_hit_counts_across_datasets_d',
    'e': 'four_strain_antibacterial_hit_overlap_venn_e',
    'f': 'growth_inhibition_log_density_violins_f',
    'g': 'TPSA_amine_count_LogP_all_strain_hits_g',
    'h': 'molecular_chemical_space_8k_coadd_tsne_h',
    'i': 'molecular_chemical_space_8k_broad39k_ecbd_tsne_i',
}
SVG_NAMES = {
    'a': 'Figure 1a (Screening workflow).svg',
    'b': 'Figure 1b (Replicate correlations).svg',
    'c': 'Figure 1c (Replicate correlation insets).svg',
    'd': 'Figure 1d (Dataset hit counts).svg',
    'e': 'Figure 1e (Hit overlap).svg',
    'f': 'Figure 1f (Growth inhibition distributions).svg',
    'g': 'Figure 1g (TPSA amines LogP).svg',
    'h': 'Figure 1h (8K COADD chemical space).svg',
    'i': 'Figure 1i (8K Broad39K ECBD chemical space).svg',
}
# PDF coordinates, measured from the top left of the current native artboard.
PANEL_BOXES = {
    'a': (0, 0, 423, 113), 'b': (0, 113, 336, 480),
    'c': (336, 113, 428, 480), 'd': (423, 0, 655, 123),
    'e': (655, 0, 875.015, 123), 'f': (423, 121, 607, 288),
    'g': (607, 121, 875.015, 288), 'h': (423, 288, 677, 539.274),
    'i': (650, 288, 875.015, 539.274),
}

METHODS_DIRECTORY = Path('/Users/asselism/Desktop/Collins_Lab/Repositories/AI-Guided_Antipseudomonal_Discovery/AI-Guided_Antipseudomonal_Discovery/figures/figure_creation_scripts')
ASSET_DIRECTORY = METHODS_DIRECTORY / 'shared/figure1_current'
INPUT_DIRECTORY = ASSET_DIRECTORY / 'inputs'
DATA_DIRECTORY = Path('/Users/asselism/Documents/Codex/2026-09-21/please-create-an-html-file-that/outputs/final_figure1_source_reproduction/1d')
OUTPUT_DIRECTORY = DATA_DIRECTORY
