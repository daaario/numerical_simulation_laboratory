#!/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
import colorsys
import sys
import glob
from pathlib import Path

# plot_against_theory.py - Plots all properties' final global averages (with error) vs T, 
#                          over their analytical theoretical graph.
# Usage: ./plot_against_theory.py <directory_containing_T_sweep_results>

if len(sys.argv) == 1:
    output_directory = 'output'
elif len(sys.argv) == 2:
    output_directory = sys.argv[1]
else:
    print("""
plot_against_theory.py - Plots all properties' final global averages (with error) vs T, 
                         over their analytical theoretical graph.
Usage: ./plot_against_theory.py <directory_containing_T_sweep_results>
""")
    sys.exit(1)

PROPERTIES = [
    ('e', 'energy_per_particle.csv', 'Energy per particle', '$e$'),
    ('c', 'specific_heat.csv', 'Specific heat', '$c_V$'),
    ('m', 'magnetization.csv', 'Magnetization', '$m$'),
    ('x', 'susceptibility.csv', 'Susceptibility', "$\\chi$")
]

# Everything colored in this figure is derived from just two anchor tints,
# so panel colors, curve gradients and the deviation scale all read as one family.
ANCHOR_COOL = '#5D7EA8'   # dusty blue   -> low end (energy panel / negative deviation)
ANCHOR_WARM = '#B4708A'   # dusty rose   -> high end (susceptibility panel / positive deviation)
NEUTRAL = '#8A8F98'

def _hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4))

def _rgb_to_hex(rgb):
    return '#' + ''.join(f'{int(round(min(max(c, 0), 1) * 255)):02x}' for c in rgb)

def _adjust_lightness(hex_color, factor):
    """factor > 1 lightens toward white, factor < 1 darkens toward black (in HLS space)."""
    h, l, s = colorsys.rgb_to_hls(*_hex_to_rgb(hex_color))
    l = l + (1 - l) * (factor - 1) if factor >= 1 else l * factor
    return _rgb_to_hex(colorsys.hls_to_rgb(h, min(max(l, 0), 1), s))

_duotone = LinearSegmentedColormap.from_list('duotone', [ANCHOR_COOL, ANCHOR_WARM])
_positions = np.linspace(0, 1, len(PROPERTIES))

# panel identity: a near-white tint (background), a light and a dark tone (curve gradient)
_NEUTRAL_SHADES = {
    'tint': _adjust_lightness(NEUTRAL, 1.9),
    'light': _adjust_lightness(NEUTRAL, 1.3),
    'dark': _adjust_lightness(NEUTRAL, 0.75),
}
PALETTE = {key: _NEUTRAL_SHADES for key, _, _, _ in PROPERTIES}

# deviation-from-theory scale (in units of sigma): same two anchors, white at zero
PULL_CMAP = LinearSegmentedColormap.from_list('pull', [ANCHOR_COOL, ANCHOR_WARM])
PULL_NORM = Normalize(vmin=0, vmax=3, clip=True)

def load_data():
    rows = []
    for dir in sorted(glob.glob(f"{output_directory}/T_*")):
        # read T value from the logfile
        logfile_path = Path(dir) / 'logfile.log'
        if not logfile_path.exists():
            sys.exit(f"Error: {logfile_path} does not exist.")
        with open(logfile_path, 'r') as logfile:
            lines = [l for l in logfile if 'TEMPERATURE:' in l]
            if not lines:
                sys.exit(f"Error: {logfile_path} does not contain any line matching 'TEMPERATURE:'.")
            T = float(lines[0].split()[1])

        # read the final global averages from each observable's file
        row = {'T': T}
        for key, filename, _, _ in PROPERTIES:
            df = pd.read_csv(Path(dir) / filename, sep=',')
            row[key] = df.iloc[-1, 3]
            row[key + '_err'] = np.maximum(df.iloc[-1, 4], 0)

        rows.append(row)
    return pd.DataFrame(rows).sort_values('T').reset_index(drop=True)

def load_params():
    params = {}
    directories = sorted(glob.glob(f"{output_directory}/T_*"))
    for key, match, logfile_index in [('Tmin', 'TEMPERATURE:', 0), ('Tmax', 'TEMPERATURE:', -1), ('J', 'J:', 0), ('H', 'H:', 0), ('N', 'NPART:', 0)]:
        with open(Path(directories[logfile_index]) / 'logfile.log', 'r') as logfile:
            lines = [l for l in logfile if match in l]
            if not lines:
                sys.exit(f"Error: no lines containing '{match}' in {logfile}")
            params[key] = float(lines[0].split()[1])

    return params

def theoretical_values(T, J, H, N):
    """Analytical finite-N, 1D Ising chain values at temperature(s) T (scalar or array)."""
    T = np.asarray(T, dtype=float)
    beta = 1 / T
    tanh = np.tanh(beta * J)
    tanhN = tanh**N
    coth = 1 / tanh
    sinh = np.sinh(beta * H)
    cosh = np.cosh(beta * H)

    # magnetization, also valid for nonzero H
    l1 = np.exp(beta * J) * cosh + np.sqrt(np.exp(2 * beta * J) * cosh * cosh - 2 * np.sinh(2 * beta * J))
    l2 = np.exp(beta * J) * cosh - np.sqrt(np.exp(2 * beta * J) * cosh * cosh - 2 * np.sinh(2 * beta * J))
    l2l1 = l2 / l1
    C = np.exp(beta * J) * cosh / np.sqrt(np.exp(2 * beta * J) * cosh * cosh - 2 * np.sinh(2 * beta * J))
    # magnetization
    m = np.exp(beta * J) * sinh * ((1 + C) / l1 + (l2l1**N) * (1 - C) / (l2 * l2)) / (1 + l2l1**N)
    # energy
    e = -J * (tanh + coth * tanhN) / (1 + tanhN) - H * m
    # susceptibility
    x = beta * np.exp(2 * beta * J) * (1 - tanhN) / (1 + tanhN)
    # specific heat
    c = ((beta * J)**2) * ((1 + tanhN + (N - 1) * tanh**2 + (N - 1) * coth**2 * tanhN) / (1 + tanhN)
                            - N * ((tanh + coth * tanhN) / (1 + tanhN))**2)

    return {'e': e, 'c': c, 'm': m, 'x': x}

def theoretical_predictions(params):
    T = np.linspace(params['Tmin'], params['Tmax'], num=200)
    values = theoretical_values(T, params['J'], params['H'], params['N'])
    return {key: (T, values[key]) for key in values}

def gradient_curve(ax, x, y, light, dark, lw=2.6):
    """Draw (x, y) as a line whose color eases from `light` to `dark` along x."""
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    cmap = LinearSegmentedColormap.from_list('grad', [light, dark])
    lc = LineCollection(segments, cmap=cmap, norm=Normalize(x.min(), x.max()),
                         linewidth=lw, capstyle='round', zorder=2)
    lc.set_array(x[:-1])
    ax.add_collection(lc)

def plot(data, functions, theory_at_data):
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'text.color': '#3A3A3A',
        'axes.labelcolor': '#3A3A3A',
        'xtick.color': '#3A3A3A',
        'ytick.color': '#3A3A3A',
    })

    fig, axs = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    fig.patch.set_facecolor('white')
    fig.suptitle('Ising 1D (chain) - Simulated values vs analytical predictions',
                 fontsize=14, fontweight='medium', color='#2E2E2E')

    scatter = None
    for ax, (key, _, title, label) in zip(axs.flat, PROPERTIES):
        colors = PALETTE[key]
        T, y = functions[key]

        ax.set_facecolor(colors['tint'])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(True, color='white', linewidth=1.3, zorder=0)
        ax.tick_params(length=0)

        gradient_curve(ax, T, y, colors['light'], colors['dark'])

        # pull: how many error-bars away each point sits from the theoretical value
        pull = np.abs( data[key].to_numpy() - theory_at_data[key]) / data[key + '_err'].to_numpy()
        ax.errorbar(data['T'], data[key], yerr=data[key + '_err'], fmt='none',
                    ecolor='#9A9A9A', elinewidth=1.0, capsize=3, zorder=3)
        scatter = ax.scatter(data['T'], data[key], c=pull, cmap=PULL_CMAP, norm=PULL_NORM,
                              edgecolor='white', linewidth=1.1, s=60, zorder=4)

        ax.set_title(title, fontsize=11, color='#2E2E2E', pad=8)
        ax.set_xlabel('$T$')
        ax.set_ylabel(label)
        ax.set_xlim(T.min(), T.max())

    fig.colorbar(scatter, ax=axs.ravel().tolist(), label='deviation from theory (σ)', shrink=0.7)

    return fig

def main():
    data = load_data()
    params = load_params()
    functions = theoretical_predictions(params)
    theory_at_data = theoretical_values(data['T'].to_numpy(), params['J'], params['H'], params['N'])
    plot(data, functions, theory_at_data)
    plt.savefig(f"{output_directory}/measures_vs_T.png")
    plt.show()

if __name__ == '__main__':
    main()
