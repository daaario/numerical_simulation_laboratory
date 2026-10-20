#!/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import colorsys
import sys
import glob
import math
from pathlib import Path

# plot_pofv.py - Plots the velocity distribution P(v) for every temperature in a
#                T-sweep: global average + error bars (bars), overlaid with the
#                per-block estimates (thin gray lines) as a fluctuation check.
#                Panels are arranged in a grid, one per T, colored by T.
# Usage: ./plot_pofv.py <directory_containing_T_sweep_results>

if len(sys.argv) ==1:
    output_directory = 'output'
elif len(sys.argv) == 2:
    output_directory = sys.argv[1]
else:
    print("""
plot_pofv.py - Plots P(v) global average (with error bars) overlaid with per-block
               fluctuations, for every T in a sweep, arranged in a grid.
Usage: ./plot_pofv.py <directory_containing_T_sweep_results>
""")
    sys.exit(1)


# --- adjust these to match your simulation's setup ---
VMIN = 0.0                     # lower edge of the histogram range used when binning v
M = 1.0                        # particle mass (LJ reduced units)
GLOBAL_FILE = 'pofv_global.csv'
BLOCKS_FILE = 'pofv_blocks.csv'
# -------------------------------------------------------------
# VMAX is set from the data itself once all runs are loaded, as
# VMAX = 5.0 * sqrt(T_max / M), so it isn't a fixed constant here.

# Same two-anchor identity as plot_ising.py: here the anchors sweep across T
# itself (cool = low T, warm = high T), so panel color doubles as a legend.
ANCHOR_COOL = '#5D7EA8'   # dusty blue   -> low T
ANCHOR_WARM = '#B4708A'   # dusty rose   -> high T


def _hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _adjust_lightness(hex_color, factor):
    """factor > 1 lightens toward white, factor < 1 darkens toward black (in HLS space)."""
    h, l, s = colorsys.rgb_to_hls(*_hex_to_rgb(hex_color))
    l = l + (1 - l) * (factor - 1) if factor >= 1 else l * factor
    return colorsys.hls_to_rgb(h, min(max(l, 0), 1), s)


def _rgba_to_hex(rgba):
    return '#%02x%02x%02x' % tuple(int(round(min(max(c, 0), 1) * 255)) for c in rgba[:3])


T_CMAP = LinearSegmentedColormap.from_list('t_scale', [ANCHOR_COOL, ANCHOR_WARM])


def load_runs():
    """Load, for each T_* directory: T, the global P(v)+error, and the per-block P(v)."""
    runs = []
    for d in sorted(glob.glob(f"{output_directory}/T_*")):
        logfile_path = Path(d) / 'logfile.log'
        if not logfile_path.exists():
            sys.exit(f"Error: {logfile_path} does not exist.")
        with open(logfile_path, 'r') as logfile:
            lines = [l for l in logfile if 'TEMPERATURE:' in l]
            if not lines:
                sys.exit(f"Error: {logfile_path} does not contain a line matching 'TEMPERATURE:'.")
            T = float(lines[0].split()[1])

        glob_path = Path(d) / GLOBAL_FILE
        blocks_path = Path(d) / BLOCKS_FILE
        if not glob_path.exists() or not blocks_path.exists():
            sys.exit(f"Error: {glob_path} or {blocks_path} does not exist.")

        glob_df = pd.read_csv(glob_path)
        blocks_df = pd.read_csv(blocks_path)
        bin_cols = [c for c in blocks_df.columns if c.startswith('BIN')]

        runs.append({'T': T, 'global': glob_df, 'blocks': blocks_df, 'bin_cols': bin_cols})

    if not runs:
        sys.exit(f"Error: no T_* directories found in {output_directory}.")

    return sorted(runs, key=lambda r: r['T'])


def grid_shape(n):
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    return rows, cols


def plot(runs):
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'text.color': '#3A3A3A',
        'axes.labelcolor': '#3A3A3A',
        'xtick.color': '#3A3A3A',
        'ytick.color': '#3A3A3A',
    })

    n = len(runs)
    rows, cols = grid_shape(n)
    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.6 * rows),
                             constrained_layout=True, squeeze=False)
    fig.patch.set_facecolor('white')
    fig.suptitle('Lennard-Jones — Velocity distribution $P(v)$ across $T$',
                 fontsize=14, fontweight='medium', color='#2E2E2E')

    Ts = np.array([r['T'] for r in runs])
    t_norm = Normalize(vmin=Ts.min(), vmax=Ts.max())
    y_max = max(r['blocks'][r['bin_cols']].to_numpy().max() for r in runs) * 1.1

    T_max = Ts.max()
    vmax = 5.0 * math.sqrt(T_max / M)

    for ax, r in zip(axs.flat, runs):
        n_bins = len(r['global'])
        edges = np.linspace(VMIN, vmax, n_bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        width = edges[1] - edges[0]

        color = T_CMAP(t_norm(r['T']))
        tint = _adjust_lightness(_rgba_to_hex(color), 1.9)

        ax.set_facecolor(tint)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(True, color='white', linewidth=1.2, zorder=0)
        ax.tick_params(length=0)

        # plot 1 (global average + error), the headline result for this T
        ax.bar(centers, r['global']['GLOBAL_AV'], width=width * 0.9,
               yerr=r['global']['ERROR'], capsize=3,
               color=color, edgecolor='white', linewidth=0.8,
               error_kw={'ecolor': '#5A5A5A', 'elinewidth': 1.0}, zorder=2)

        # plot 2 (per-block fluctuations), drawn last so plot 1 sits on top
        for _, row in r['blocks'].iterrows():
            ax.plot(centers, row[r['bin_cols']].to_numpy(dtype=float),
                     color='#9A9A9A', alpha=0.35, lw=0.9, zorder=1)

        ax.set_title(f"$T = {r['T']:.3g}$", fontsize=11, color='#2E2E2E', pad=6)
        ax.set_xlabel('$v$')
        ax.set_ylabel('$P(v)$')
        ax.set_xlim(VMIN, vmax)
        ax.set_ylim(0, y_max)

    # hide any unused panels if the grid isn't fully filled (e.g. 7 runs -> 3x3 grid)
    for ax in axs.flat[n:]:
        ax.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=T_CMAP, norm=t_norm)
    fig.colorbar(sm, ax=axs.ravel().tolist(), label='$T$', shrink=0.7)

    return fig


def main():
    runs = load_runs()
    plot(runs)
    plt.savefig(f"{output_directory}/pofv_vs_T.png", dpi=150)
    plt.show()


if __name__ == '__main__':
    main()
