from pathlib import Path

ANALYSIS_DIRECTORY = Path(__file__).resolve().parent
INPUT_DIRECTORY = ANALYSIS_DIRECTORY / "inputs"
INPUT_PATHS = {"OD600": "", "Resazurin_30min": "", "Resazurin_60min": ""}
INPUT_FILENAMES = {
    "OD600": "Princeton_Bio_WTAtoD_d6EtoH.txt",
    "Resazurin_30min": "Princeton_Bio_WTAtoD_d6EtoH_resazurin_30mins.txt",
    "Resazurin_60min": "Princeton_Bio_WTAtoD_d6EtoH_resazurin_60mins.txt",
}
HIGHEST_CONCENTRATION_UM = 255.0
DILUTION_FACTOR = 2.0
COLORS = {"WT": "#00a650", "Delta6": "#a860a5"}
MARKERS = {"WT": "o", "Delta6": "s"}
FONT_PATH = "/System/Library/Fonts/Helvetica.ttc"
PNG_DPI = 300
