"""Inputs, colors, and export settings for Figure 2."""

from pathlib import Path

METHODS_DIRECTORY = Path(__file__).resolve().parent
DATA_MODE = "saved"
DATASET_PATHS = {"COADD": "", "8K": ""}
ASSET_PATH = METHODS_DIRECTORY / "shared" / "assets" / "figure2" / "figure_2_vector_layout.pdf"
OUTPUT_DIRECTORY = METHODS_DIRECTORY / "full_figure_2"
PNG_DPI = 180
WORKER_THREADS = 8
INFERENCE_BATCH_SIZE = 4096
PANEL_NAMES = {
    "a": "deepflux_gi60_gi80_delta5_hits_vs_inactive_roc_a",
    "b": "ranked_compound_delta6_and_wild_type_hit_retrieval_b",
    "c": "compound_ic50_in_wild_type_and_delta6_strains_c",
    "d": "gaussian_property_distributions_delta6_only_and_wild_type_hits_d",
    "e": "deepflux_efflux_score_compound_predictions_and_wild_type_hits_e",
    "f": "molecular_property_distributions_by_growth_inhibition_f",
    "_unlabeled_growth_inhibition": "growth_inhibition_cumulative_distributions_figure_2_unlabeled",
    "_unlabeled_compound_examples": "compound_structures_and_dose_responses_figure_2_unlabeled",
}
PANEL_BOXES = {
    "a": (0, 0, 158, 173),
    "b": (157, 0, 367, 174),
    "c": (359, 0, 657, 192),
    "d": (651, 0, 841.89, 197),
    "e": (0, 172, 479, 408),
    "f": (481, 197, 841.89, 424),
    "_unlabeled_growth_inhibition": (43, 422, 257, 595.276),
    "_unlabeled_compound_examples": (279, 424, 832, 595.276),
}
FIGURE_NAME = "deepflux_predictions_and_hit_properties_figure_2"
