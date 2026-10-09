"""Panel 1g box-and-whisker plots with the current descriptors and display limits."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/matplotlib_cache')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from violin_panels import PROPERTY_LIMITS, PROPERTY_TICKS


def plot_property_boxplots(p, props, groups, colors, save, tables):
    names = ['Heavy atoms', 'logP', 'TPSA (Å²)', 'log$_{10}$ IPC', 'Number of amines', 'Aromatic rings']
    fig = plt.figure(figsize=(257 / 72, 164 / 72))
    counts = [int(p.Category.eq(g).sum()) for g in groups]
    fig.text(.01, .96, 'Physicochemical properties', fontsize=8, weight='bold', va='top')
    for x, label, n, col in zip([.01, .30, .67], ['WT hits', 'Δ6-only hits', 'Inactive'], counts, colors):
        fig.text(x, .86, '●', color=col, fontsize=7)
        fig.text(x + .04, .86, f'{label} ({n:,})', fontsize=6)
    rng = np.random.default_rng(42)
    audit = []
    for i, (prop, name) in enumerate(zip(props, names)):
        ax = fig.add_axes([.10 + (i % 3) * .33, .51 - (i // 3) * .40, .22, .22])
        lo, hi = PROPERTY_LIMITS[prop]
        for j, (group, col) in enumerate(zip(groups, colors)):
            a = p.loc[p.Category.eq(group), prop].to_numpy()
            assert np.isfinite(a).all()
            # Quartiles and whiskers use all observations before applying the display zoom.
            ax.boxplot([a], positions=[j], widths=.52, patch_artist=True,
                       whis=1.5, showfliers=True, manage_ticks=False,
                       boxprops={'linewidth': .5, 'facecolor': col, 'edgecolor': '#444444', 'alpha': .55},
                       medianprops={'color': 'black', 'linewidth': .75},
                       whiskerprops={'linewidth': .5, 'color': '#444444'},
                       capprops={'linewidth': .5, 'color': '#444444'},
                       flierprops={'marker': 'o', 'markersize': 1.3, 'markerfacecolor': col,
                                   'markeredgewidth': 0, 'alpha': .35})
            if j == 0:
                ax.scatter(j + rng.uniform(-.17, .17, len(a)), a, s=3,
                           color=col, edgecolors='black', linewidths=.15, zorder=5)
            below, above = int((a < lo).sum()), int((a > hi).sum())
            if below:
                ax.plot(j, lo, 'v', markersize=2.6, color=col, clip_on=False, zorder=6)
            if above:
                ax.plot(j, hi, '^', markersize=2.6, color=col, clip_on=False, zorder=6)
            q1, median, q3 = np.quantile(a, [.25, .5, .75])
            inside = a[(a >= q1 - 1.5 * (q3 - q1)) & (a <= q3 + 1.5 * (q3 - q1))]
            audit.append({'Property': prop, 'Group': group, 'N': len(a),
                          'Q1': q1, 'Median': median, 'Q3': q3,
                          'Lower_whisker': inside.min(), 'Upper_whisker': inside.max(),
                          'Display_min': lo, 'Display_max': hi,
                          'Below_display': below, 'Above_display': above,
                          'Percent_in_display': 100 * ((a >= lo) & (a <= hi)).mean(),
                          'Actual_min': a.min(), 'Actual_max': a.max()})
        ax.set_xlim(-.65, 2.65)
        ax.set_ylim(lo, hi)
        ax.set_title(name, fontsize=6.5, pad=3)
        ax.set_xticks([0, 1, 2], ['WT', 'Δ6', 'I'])
        ax.tick_params(axis='x', length=0)
        ax.set_yticks(PROPERTY_TICKS[prop])
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=6, pad=1.5)
    fig.text(.01, .035, 'Δ6 = Δ6-only; I = inactive in WT and Δ6; ▲/▼ beyond axis', fontsize=5.5)
    save(fig, 'g_physicochemical_properties')
    pd.DataFrame(audit).to_csv(tables / 'property_boxplot_display_audit.csv', index=False)
    return counts
