"""Score CO-ADD with both supplied ten-member DeepFlux ensembles."""

from pathlib import Path
import sys

METHODS_DIRECTORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(METHODS_DIRECTORY))
from shared.figure2_analysis import score_coadd

if __name__ == "__main__":
    score_coadd()
