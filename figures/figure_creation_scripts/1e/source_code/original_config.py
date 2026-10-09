"""Input locations, analysis settings, and output locations for Figure 1."""

from pathlib import Path

METHODS_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_DIRECTORY = METHODS_DIRECTORY / "full_figure"
CACHE_DIRECTORY = METHODS_DIRECTORY / "shared" / "computed_data"
ASSET_DIRECTORY = METHODS_DIRECTORY / "shared" / "assets"
SAVED_DATA_DIRECTORY = METHODS_DIRECTORY / "shared" / "saved_plot_data"

DATA_MODE = "saved"
RECOMPUTE_EMBEDDINGS = False
DATASET_PATHS = {
    "8K": "",
    "39K": "",
    "ECBD": "",
    "COADD": "",
}

STRAINS = ("WT", "Delta40", "Delta3", "Delta6")
RANDOM_SEED = 42
WORKER_THREADS = 8
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048
PCA_COMPONENTS = 50
TSNE_PERPLEXITY = 30
PNG_DPI = 180

PANEL_NAMES = {
    "a": "antipseudomonal_compound_screening_workflow_a",
    "b": "replicate_growth_inhibition_correlations_across_four_strains_b",
    "c": "strain_specific_antibacterial_hit_counts_c",
    "d": "strain_specific_antibacterial_hit_intersections_upset_d",
    "e": "four_strain_antibacterial_hit_overlap_venn_e",
    "f": "molecular_chemical_space_8k_broad39k_ecbd_tsne_f",
    "g": "molecular_chemical_space_8k_coadd_tsne_g",
}
FIGURE_NAME = "antipseudomonal_screening_and_chemical_space_figure_1"

PANEL_DIRECTORIES = {panel: METHODS_DIRECTORY / f"1{panel}" for panel in PANEL_NAMES}

PANEL_BOXES = {
    "a": (0, 350, 385, 469.274),
    "b": (0, 0, 350, 350),
    "c": (392, 334, 520, 469.274),
    "d": (522, 334, 657, 469.274),
    "e": (658, 340, 870, 469.274),
    "f": (361, 28, 622, 332),
    "g": (625, 28, 868, 332),
}
COLORS_CMYK = {
    "WT": (1, 0, 1, 0),
    "Delta40": (0, 0.418188750743866, 1, 0),
    "Delta3": (0.177340358495712, 0.909941256046295, 0.953368425369263, 0.070771344006062),
    "Delta6": (0.360784313725490, 0.737254901960784, 0.011764705882353, 0),
    "Delta5": (0.578125, 0.04296875, 0.08203125, 0),
}
