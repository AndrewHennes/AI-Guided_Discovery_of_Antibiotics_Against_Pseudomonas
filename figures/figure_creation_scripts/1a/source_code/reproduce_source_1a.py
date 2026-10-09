"""Reproduce the source artwork and calculations associated with Figure 1a."""
from pathlib import Path
import sys

SOURCE_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE_DIRECTORY))
import config
from shared.rendering import render_panel

if __name__ == '__main__':
    config.OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    render_panel(config.SOURCE_PANEL)
