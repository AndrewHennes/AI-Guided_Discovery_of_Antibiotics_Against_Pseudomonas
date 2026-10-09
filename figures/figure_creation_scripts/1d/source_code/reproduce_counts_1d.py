"""Recalculate the four-dataset hit counts and draw the Figure 1d plot."""
from pathlib import Path
import sys

SOURCE_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE_DIRECTORY))
from shared.figure1_current_analysis import counts, plt

if __name__ == '__main__':
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': .6,
                         'xtick.major.width': .6, 'ytick.major.width': .6,
                         'xtick.major.size': 2, 'ytick.major.size': 2,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    counts()
