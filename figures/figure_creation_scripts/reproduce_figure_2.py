"""Reproduce Figure 2 and its individual panel exports."""

from pathlib import Path
import sys

METHODS_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(METHODS_DIRECTORY))
from shared.figure2_rendering import render_figure

if __name__ == "__main__":
    render_figure()
