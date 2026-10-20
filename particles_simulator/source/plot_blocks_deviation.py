#!/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import colorsys
import sys
import glob
from pathlib import Path

# plot_blocks_deviation.py - Plots all properties' block average for all computed 
#                            values of T, centered around their theoretical value.
# Usage: ./plot_blocks_deviation.py <directory_containing_T_sweep_results>

if len(sys.argv) == 1:
    output_directory = 'output'
elif len(sys.argv) == 2:
    output_directory = sys.argv[1]
else:
    print("""
plot_blocks_deviation.py - Plots all properties' block average for all computed 
                           values of T, centered around their theoretical value.
Usage: ./plot_blocks_deviation.py <directory_containing_T_sweep_results>
""")
    sys.exit(1)

PROPERTIES = [
    ('e', 'energy_per_particle.csv', 'Energy per particle'),
    ('c', 'specific_heat.csv', 'Specific heat'),
    ('m', 'magnetization.csv', 'Magnetization'),
    ('x', 'susceptibility.csv', 'Susceptibility'),
]

# same two anchors as plot.py's residual scale - here they're doing double duty
# on purpose, since temperature genuinely has a cold/hot association.
ANCHOR_COOL = '#5D7EA8'   # dusty blue -> low temperature
ANCHOR_WARM = '#B4708A'   # dusty rose -> high temperature
NEUTRAL = '#8A8F98'       # soft slate gray, shared panel background

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

TEMP_CMAP = LinearSegmentedColormap.from_list('temperature', [ANCHOR_COOL, ANCHOR_WARM])
_NEUTRAL_TINT = _adjust_lightness(NEUTRAL, 1.9)
_NEUTRAL_DARK = _adjust_lightness(NEUTRAL, 0.6)
MATHLABELS = {'e': r'$\Delta e$', 'c': r'$\Delta c$', 'm': r'$\Delta m$', 'x': r'$\Delta \chi$'}

def read_parameters(logfile_path):
    with open(logfile_path, 'r') as logfile:
        all_lines = logfile.readlines()
    params = {}
    for key, match in [('T', 'TEMPERATURE:'), ('J', 'J:'), ('H', 'H:'), ('N', 'NPART:')]:
        lines = [l for l in all_lines if match in l]
        if not lines:
            sys.exit(f"Error: {logfile_path} does not contain any line matching '{match}'.")
        params[key] = float(lines[0].split()[1])
    return params

def load_data():
    data = []

    for dir in sorted(glob.glob(f"{output_directory}/T_*")):
        logfile_path = Path(dir) / 'logfile.log'
        if not logfile_path.exists():
            sys.exit(f"Error: {logfile_path} does not exist.")
        params = read_parameters(logfile_path)
        # params is {'T': *, 'J': 1.0, 'H': 0.02, 'N': 200}

        series = {}
        for key, filename, _ in PROPERTIES:
            csv_path = Path(dir) / filename
            if not csv_path.exists():
                sys.exit(f"Error: {csv_path} does not exist.")
            df = pd.read_csv(csv_path, sep=',')
            series[key] = (df.iloc[:, 0].to_numpy(), df.iloc[:, 1].to_numpy())
        # series is {'e': ([1, 2, 3, ...], (39.02, 39.17, 38.98, ...), 'c': ([...], [...]), ...}
        
        data.append({'T': params['T'], 'J': params['J'], 'H': params['H'], 'N': params['N'], 'series': series})

    if not data:
        sys.exit(f"Error: no '{output_directory}/T_*' directories found.")

    return data

def theoretical_value(prop, T, J, H, N):
    beta = 1/T
    tanh = np.tanh(beta * J)

    if prop == 'e':
        tanhN = tanh**N
        coth = 1 / tanh
        sinh = np.sinh(beta * H)
        cosh = np.cosh(beta * H)
        l1 = np.exp(beta * J)*cosh + np.sqrt( np.exp(2*beta*J)*cosh*cosh - 2*np.sinh(2*beta*J) )
        l2 = np.exp(beta * J)*cosh - np.sqrt( np.exp(2*beta*J)*cosh*cosh - 2*np.sinh(2*beta*J) )
        l2l1 = l2/l1
        C = np.exp(beta * J)*cosh / np.sqrt(np.exp(2*beta*J)*cosh*cosh - 2*np.sinh(2*beta*J))
        return -J * (tanh + coth*tanhN) / (1 + tanhN) - H *( np.exp(beta * J) * sinh * ( (1 + C)/l1 + (l2l1**N)*(1 - C)/(l2*l2) ) / ( 1 + l2l1**N ) )
    if prop == 'c':
        tanhN = tanh**N
        coth = 1 / tanh
        return ((beta * J)**2) * ( (1 + tanhN + (N-1)*(tanh**2) + (N-1)*(coth**2)*tanhN)/(1+tanhN) - N*((tanh + coth*tanhN)/(1 + tanhN))**2)
    if prop == 'm':
        sinh = np.sinh(beta * H)
        cosh = np.cosh(beta * H)
        l1 = np.exp(beta * J)*cosh + np.sqrt( np.exp(2*beta*J)*cosh*cosh - 2*np.sinh(2*beta*J) )
        l2 = np.exp(beta * J)*cosh - np.sqrt( np.exp(2*beta*J)*cosh*cosh - 2*np.sinh(2*beta*J) )
        l2l1 = l2/l1
        C = np.exp(beta * J)*cosh / np.sqrt(np.exp(2*beta*J)*cosh*cosh - 2*np.sinh(2*beta*J))
        return np.exp(beta * J) * sinh * ( (1 + C)/l1 + (l2l1**N)*(1 - C)/(l2*l2) ) / ( 1 + l2l1**N )
    if prop == 'x':
        tanhN = tanh**N
        return beta * np.exp(2 * beta * J) * (1 - tanhN) / (1 + tanhN)

    raise ValueError(prop)

def plot_blocks(data):
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'text.color': '#3A3A3A',
        'axes.labelcolor': '#3A3A3A',
        'xtick.color': '#3A3A3A',
        'ytick.color': '#3A3A3A',
    })

    T = np.array([d['T'] for d in data])
    norm = Normalize(vmin=T.min(), vmax=T.max())

    fig, axs = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    fig.patch.set_facecolor('white')
    fig.suptitle('Ising 1D (chain) - block average deviation from theory',
                 fontsize=14, fontweight='medium', color='#2E2E2E')

    for ax, (key, _, title) in zip(axs.flat, PROPERTIES):
        ax.set_facecolor(_NEUTRAL_TINT)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(True, color='white', linewidth=1.3, zorder=0)
        ax.tick_params(length=0)
        ax.axhline(0, color=_NEUTRAL_DARK, linewidth=1.2, linestyle='--', alpha=0.6, zorder=1)

        for run in sorted(data, key=lambda r: r['T']):  # low-T lines drawn first, so higher-T lines consistently layer on top
            blocks, block_avgs = run['series'][key]
            theory = theoretical_value(key, run['T'], run['J'], run['H'], run['N'])
            deviation = block_avgs - theory  # new array - don't mutate the series stored in `data`

            ax.plot(blocks, deviation, color=TEMP_CMAP(norm(run['T'])), linewidth=1.8, alpha=0.9, zorder=2)

        ax.set_title(title, fontsize=11, color='#2E2E2E', pad=8)
        ax.set_xlabel('block')
        ax.set_ylabel(MATHLABELS[key])

    sm = plt.cm.ScalarMappable(cmap=TEMP_CMAP, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=axs.ravel().tolist(), label='temperature ($T$)', shrink=0.8)

    return fig

def main():
    data = load_data()
    plot_blocks(data)
    plt.savefig(f"{output_directory}/block_deviations.png", format='png', dpi=900)
    plt.show()

if __name__ == '__main__':
    main()
