"""Violin panels for Figure 1, with explicit density and display scales."""
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/matplotlib_cache')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

SCREEN_N = 7994
DENSITY_FLOOR = 1e-5  # Probability density per percentage point.
INHIBITION_BANDWIDTH = 2.0
PROPERTY_LIMITS = {
    'MW': (150, 550), 'logP': (0, 6.5), 'TPSA': (0, 150),
    'Log10_IPC': (1.5, 9.5), 'Rotatable_bonds': (-.5, 10.5),
    'Aromatic_rings': (-.4, 5.4), 'Amine_count': (-.4, 4.4), 'Heavy_atoms': (5, 40), 'Alcohol_count': (-.3, 3.3),
}
PROPERTY_TICKS = {
    'MW': [200, 350, 500], 'logP': [0, 3, 6], 'TPSA': [0, 75, 150],
    'Log10_IPC': [2, 5, 8], 'Rotatable_bonds': [0, 5, 10],
    'Aromatic_rings': [0, 2, 4], 'Amine_count': [0, 2, 4], 'Heavy_atoms': [10, 25, 40], 'Alcohol_count': [0, 1, 2, 3],
}


def _style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=6, pad=1.5)


def plot_inhibition_violins(values, summary, strains, labels, colors, save, tables):
    x = np.linspace(-60, 110, 1361)
    densities = []
    for strain in strains:
        v = values[f'Average_{strain}_Percent_Growth_Inhibition'].dropna().to_numpy()
        kde = stats.gaussian_kde(v, bw_method=INHIBITION_BANDWIDTH / np.std(v, ddof=1))
        densities.append(kde(x))
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


def plot_property_violins(p, props, groups, colors, save, tables):
    names = ['MW (Da)', 'logP', 'TPSA (Å²)', 'log$_{10}$ IPC', 'Rotatable bonds', 'Aromatic rings']
    fig = plt.figure(figsize=(257 / 72, 164 / 72))
    counts = [int(p.Category.eq(g).sum()) for g in groups]
    fig.text(.01, .96, 'Physicochemical properties', fontsize=8, weight='bold', va='top')
    for x, g, n, col in zip([.01, .30, .67], ['WT hits', 'Δ6-only hits', 'Inactive'], counts, colors):
        fig.text(x, .86, '●', color=col, fontsize=7)
        fig.text(x + .04, .86, f'{g} ({n:,})', fontsize=6)
    rng = np.random.default_rng(42)
    audit = []
    for i, (prop, name) in enumerate(zip(props, names)):
        ax = fig.add_axes([.10 + (i % 3) * .33, .51 - (i // 3) * .40, .22, .22])
        pooled = p[prop].to_numpy()
        assert np.isfinite(pooled).all(), prop
        robust_sd = min(np.std(pooled, ddof=1), (np.quantile(pooled, .75) - np.quantile(pooled, .25)) / 1.34)
        bandwidth = .9 * robust_sd * len(pooled) ** (-.2)
        if prop in ['Rotatable_bonds', 'Aromatic_rings']:
            bandwidth = max(bandwidth, .35)
        assert bandwidth > 0
        lo, hi = PROPERTY_LIMITS[prop]
        # Evaluate the KDE over the full data support; display zoom does not filter it.
        y = np.linspace(min(pooled.min(), lo), max(pooled.max(), hi), 1201)
        for j, (group, col) in enumerate(zip(groups, colors)):
            a = p.loc[p.Category.eq(group), prop].to_numpy()
            density = stats.gaussian_kde(a, bw_method=bandwidth / np.std(a, ddof=1))(y)
            width = .39 * density / density.max()
            # Cut at observed group extrema to avoid implying observations beyond them.
            keep = (y >= a.min()) & (y <= a.max())
            yy = np.r_[a.min(), y[keep], a.max()]
            ww = np.r_[0, width[keep], 0]
            ax.fill_betweenx(yy, j - ww, j + ww, facecolor=col, alpha=.58, edgecolor=col, linewidth=.5)
            q = np.quantile(a, [.25, .5, .75])
            ax.plot([j, j], q[[0, 2]], color='#333333', linewidth=1.0, zorder=3)
            ax.plot(j, q[1], 'o', markersize=1.7, markerfacecolor='white', markeredgecolor='#333333', markeredgewidth=.4, zorder=4)
            if j == 0:
                ax.scatter(j + rng.uniform(-.17, .17, len(a)), a, s=3, color=col, edgecolors='black', linewidths=.15, zorder=5)
            below, above = int((a < lo).sum()), int((a > hi).sum())
            if below:
                ax.plot(j, lo, 'v', markersize=2.6, color=col, clip_on=False, zorder=6)
            if above:
                ax.plot(j, hi, '^', markersize=2.6, color=col, clip_on=False, zorder=6)
            audit.append({'Property': prop, 'Group': group, 'N': len(a), 'Bandwidth': bandwidth,
                          'Display_min': lo, 'Display_max': hi, 'Below_display': below, 'Above_display': above,
                          'Percent_in_display': 100 * ((a >= lo) & (a <= hi)).mean(),
                          'Actual_min': a.min(), 'Actual_max': a.max()})
        ax.set_xlim(-.65, 2.65)
        ax.set_ylim(lo, hi)
        ax.set_title(name, fontsize=6.5, pad=3)
        ax.set_xticks([0, 1, 2], ['WT', 'Δ6', 'I'])
        ax.tick_params(axis='x', length=0)
        ax.set_yticks(PROPERTY_TICKS[prop])
        _style(ax)
    fig.text(.01, .035, 'Δ6 = Δ6-only; I = inactive in WT and Δ6; ▲/▼ beyond axis', fontsize=5.5)
    save(fig, 'g_physicochemical_properties')
    pd.DataFrame(audit).to_csv(tables / 'property_violin_display_audit.csv', index=False)
    return counts
