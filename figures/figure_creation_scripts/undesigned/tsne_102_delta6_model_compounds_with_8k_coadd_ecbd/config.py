"""Input and analysis settings for the 102-compound, 8K, CO-ADD, and ECBD t-SNE."""

from pathlib import Path

DIRECTORY = Path(__file__).resolve().parent
DATASET_PATHS = {"8K": "", "COADD": "", "ECBD": "", "New102": ""}
RECOMPUTE_EMBEDDING = False
RANDOM_SEED = 42
WORKER_THREADS = 8
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048
PCA_COMPONENTS = 50
TSNE_PERPLEXITY = 30
EARLY_EXAGGERATION_ITERATIONS = 250
TSNE_ITERATIONS = 500
GROWTH_INHIBITION_CUTOFF = 80
FIGURE_STEM = "standardized_102_delta6_model_compounds_with_8k_coadd_ecbd_tsne"

PLOT_STYLE = {
    "canvas_width_pt": 243.0,
    "canvas_height_pt": 304.0,
    "plot_bounds_pt": [[6.0, 231.0], [21.0, 241.0]],
    "background_radius_pt": 0.303108739853,
    "foreground_radius_pt": 1.705688905716,
    "foreground_outline_width_pt": 0.25,
    "background_color": "#c4c4c4",
    "background_opacity": 0.5,
    "foreground_opacity": 0.949997,
    "wt_color": "#00a650",
    "delta6_color": "#a860a5",
    "outline_color": "#030405",
    "text_color": "#221f1f",
    "font_family": "Helvetica",
    "legend_font_size_pt": 8.4,
}
