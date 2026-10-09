"""Violin panels for Figure 1, with explicit density and display scales."""
import json
import os
from figure1_current_config import DATA_DIRECTORY
os.environ.setdefault('MPLCONFIGDIR', str(DATA_DIRECTORY / 'matplotlib_cache'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCREEN_N = 7994
DENSITY_FLOOR = 1e-5  # Probability density per percentage point.
INHIBITION_BANDWIDTH = 2.0
def _style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=6, pad=1.5)


def plot_inhibition_violins(values, summary, strains, labels, colors, save, tables):
    x = np.linspace(-60, 110, 1361)
    densities = []
    for strain in strains:
        v = values[f'Average_{strain}_Percent_Growth_Inhibition'].dropna().to_numpy()
        # Evaluate the fixed-bandwidth Gaussian KDE directly. This avoids the
        # platform BLAS covariance warnings from rescaling scipy.gaussian_kde.
        density = np.empty_like(x)
        for start in range(0, len(x), 128):
            z = (x[start:start + 128, None] - v[None, :]) / INHIBITION_BANDWIDTH
            density[start:start + 128] = np.exp(-.5 * z * z).mean(axis=1) / (INHIBITION_BANDWIDTH * np.sqrt(2 * np.pi))
        assert np.isfinite(density).all() and (density >= 0).all()
        densities.append(density)
    maximum = max(d.max() for d in densities)
    # A shared reference floor makes every half-width nonnegative and comparable.
    # This transforms density, not inhibition, and does not renormalize each strain.
    denominator = np.log10(maximum / DENSITY_FLOOR)
    fig = plt.figure(figsize=(185 / 72, 164 / 72))
    ax = fig.add_axes([.27, .235, .59, .56])
    density_table = {'Mean_growth_inhibition_percent': x}
    for i, (strain, density) in enumerate(zip(strains, densities)):
        base = 3 - i
        width = .405 * np.maximum(np.log10(np.maximum(density, DENSITY_FLOOR) / DENSITY_FLOOR), 0) / denominator
        ax.fill_between(x, base - width, base + width, facecolor=colors[strain], alpha=.55, linewidth=0)
        ax.plot(x, base + width, color=colors[strain], linewidth=.65)
        ax.plot(x, base - width, color=colors[strain], linewidth=.65)
        v = values[f'Average_{strain}_Percent_Growth_Inhibition'].dropna()
        q = v.quantile([.25, .5, .75]).to_numpy()
        ax.plot(q[[0, 2]], [base, base], color='black', linewidth=1.4, solid_capstyle='butt')
        ax.plot(q[1], base, 'o', markersize=2, markerfacecolor='white', markeredgecolor='black', markeredgewidth=.4)
        hits = int(summary.loc[summary.Strain.eq(strain), 'Binary_hit_structures'].iloc[0])
        ax.text(114, base, str(hits), ha='left', va='center', fontsize=6)
        density_table[f'{strain}_density'] = density
        density_table[f'{strain}_half_width'] = width
    ax.axvline(80, color='#555555', linestyle=(0, (2, 2)), linewidth=.6, zorder=0)
    ax.set_xlim(-60, 110)
    ax.set_ylim(-.48, 3.48)
    ax.set_xticks([-40, 0, 40, 80])
    ax.set_yticks([3, 2, 1, 0], [f'{labels[s]}\nScreen n = {SCREEN_N:,}' for s in strains])
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('Mean growth inhibition (%)', fontsize=7, labelpad=3)
    _style(ax)
    ax.tick_params(axis='y', labelsize=5.5)
    ax.spines['left'].set_visible(False)
    fig.text(.03, .96, 'Growth-inhibition distribution', fontsize=8, weight='bold', va='top')
    fig.text(.03, .875, 'Violin widths use a shared log$_{10}$ density scale', fontsize=5.7)
    fig.text(.87, .81, '≥80%', fontsize=6, ha='left')
    fig.text(.03, .035, 'Density floor $10^{-5}$ per percentage point', fontsize=5.7)
    save(fig, 'f_inhibition_log_density_violins')
    pd.DataFrame(density_table).to_csv(tables / 'inhibition_violin_densities.csv', index=False)
    (tables / 'inhibition_violin_settings.json').write_text(json.dumps({
        'screen_n_annotation': SCREEN_N,
        'screen_n_source': 'Original replicate-panel screen annotation.',
        'measured_unique_structures': dict(zip(summary.Strain, summary.Measured_unique_structures.astype(int))),
        'density_units': 'Probability per percentage point of growth inhibition',
        'gaussian_bandwidth_percentage_points': INHIBITION_BANDWIDTH,
        'density_floor': DENSITY_FLOOR, 'shared_peak_density': float(maximum),
        'half_width_formula': '0.405 * max(log10(density / 1e-5), 0) / log10(shared_peak_density / 1e-5)',
        'displayed_inhibition_range': [-60, 110], 'independent_strain_normalization': False,
    }, indent=2))

