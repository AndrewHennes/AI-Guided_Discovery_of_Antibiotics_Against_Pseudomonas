"""Export Figure 2b with the shared strain palette."""

from pathlib import Path
import sys

METHODS_DIRECTORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(METHODS_DIRECTORY))

from shared.figure2_rendering import render_panel


if __name__ == "__main__":
    render_panel("b")
