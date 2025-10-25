# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Visualization."""
from importlib.resources import files

MY_MPL_STYLE = files("feupy.visualization.styles") / "mystyle.mplstyle"

# from .setting import (
#     set_leg_style, 
#     set_leg_style_models, 
#     set_leg_style_datasets, 
#     setting_leg_style,
#     generate_legend_set,
#     get_kwargs_datasets,
#     get_kwargs_models,
# )

# from .spectral_energy_distribution import show_SED
# from .sky_map import show_sky_map, create_sky_map, show_ROI_sky_map
# from .counts import show_hist_counts, show_sensitivity_curve
# from .styles.build_plotting_config import generate_marker_set
# from .styles.legend import generate_markers, plot_colortable

# List of all functions and variables to be made available when importing this module
__all__ = [
#     "generate_markers",
#     "plot_colortable",
#     "set_leg_style",
#     "set_leg_style_models",
#     "set_leg_style_datasets",
#     "setting_leg_style",
#     "generate_legend_set",
#     "get_kwargs_datasets",
#     "get_kwargs_models",
#     "show_SED",
#     "show_sky_map",
#     "show_ROI_sky_map",
#     "create_sky_map",
#     "show_hist_counts",
#     "show_sensitivity_curve",
#     'generate_marker_set',
#     'set_axes',
]

# Default color palettes for visualizations

PALETTE_DEFAULT = [
    ['seagreen', 'seagreen'],
    ['blue', 'blue'],
    ['aquamarine', 'aquamarine'],
    ['salmon', 'salmon'],
    ['crimson', 'crimson'],
    ['tan', 'tan'],
    ['turquoise', 'turquoise'],
    ['indigo', 'indigo'],
    ['pink', 'pink'],
    ['chocolate', 'chocolate'],
    ['fuchsia', 'fuchsia'],
    ['green', 'green'],
    ['purple', 'purple'],
    ['aqua', 'aqua'],
    ['olive', 'olive'],
    ['sienna', 'sienna'],
    ['orchid', 'orchid'],
    ['orangered', 'orangered'],
    ['orange', 'orange'],
    ['cadetblue', 'cadetblue'],
    ['wheat', 'wheat'],
    ['khaki', 'khaki'],
    ['peru', 'peru'],
    ['teal', 'teal'],
    ['springgreen', 'springgreen'],
    ['red', 'red'],
    ['darkblue', 'darkblue'],
    ['brown', 'brown'],
    ['tomato', 'tomato'],
    ['coral', 'coral'],
    ['yellowgreen', 'yellowgreen'],
    ['skyblue', 'skyblue'],
    ['chartreuse', 'chartreuse'],
    ['plum', 'plum'],
    ['maroon', 'maroon'],
    ['silver', 'silver'],
]

# Standard matplotlib colours
PALETTE_TABLEAU = [
    ['tab:blue', 'blue'],
    ['tab:orange', 'orange'],
    ['tab:green', 'green'],
    ['tab:purple', 'purple'],
    ['tab:brown', 'brown'],
    ['tab:pink', 'pink'],
    ['tab:gray', 'grey'],
    ['tab:olive', 'olive'],
    ['tab:cyan', 'cyan'],
]

# Colorblind-friendly palettes (IBM and Wong)
PALETTE_IBM = [
    ['#648FFF', 'blue'],
    ['#785EF0', 'dark lavender'],
    ['#DC267F', 'deep pink'],
    ['#FE6100', 'dark orange'],
    ['#FFB000', 'orange'],
    ['#FDDF49', 'yellow'],
    ['#8FB327', 'lime green'],
    ['#902499', 'purple'],
    ['#17E474', 'lime'],
    ['#74CE22', 'green'],
    ['#449C8C', 'teal'],
    ['#EA83B7', 'pink'],
]

PALETTE_WONG = [
    ['#DF0A10', 'red'],
    ['#E69F00', 'orange'],
    ['#56B4E9', 'sky blue'],
    ['#009E73', 'sea green'],
    ['#F0E442', 'yellow'],
    ['#0072B2', 'blue'],
    ['#D55E00', 'dark orange'],
    ['#CC79A7', 'pink'],
    ['#80C774', 'lime green'],
    ['#A21DCA', 'purple'],
    ['#9BFB5F', 'lawn green'],
    ['#BDD1D2', 'grey'],
]

# Default markers for visualizations
MARKERS_DEFAULT = [
    ['o', 0.95, 'circle'],
    ['v', 1, 'down-pointing triangle'],
    ['^', 1, 'up-pointing triangle'],
    ['<', 1, 'left-pointing triangle'],
    ['>', 1, 'right-pointing triangle'],
    ['8', 1, 'octagon'],
    ['s', 0.9, 'square'],
    ['p', 1, 'pentagon'],
    ['P', 1.2, 'plus'],
    ['*', 1.5, 'star'],
    ['h', 1, 'hexagon'],
    ['H', 1, 'rotated hexagon'],
    ['X', 1, 'cross'],
    ['D', 0.75, 'diamond'],
    ['d', 0.9, 'thin diamond'],
]

MARKERS_DEFAULT_DICT = {
    'o': [0.95, 'circle'],
    'v': [1, 'down-pointing triangle'],
    '^': [1, 'up-pointing triangle'],
    '<': [1, 'left-pointing triangle'],
    '>': [1, 'right-pointing triangle'],
    '8': [1, 'octagon'],
    's': [0.9, 'square'],
    'p': [1, 'pentagon'],
    'P': [1.2, 'plus'],
    '*': [1.5, 'star'],
    'h': [1, 'hexagon'],
    'H': [1, 'rotated hexagon'],
    'X': [1, 'cross'],
    'D': [0.75, 'diamond'],
    'd': [0.9, 'thin diamond']
}
# MARKERS = ['H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h']

# Default line styles for visualizations
LINESTYLES_DEFAULT = [
    'solid', 
    (0, (1, 10)),  # loosely dotted
    'dotted',
    (0, (1, 1)),  # densely dotted
    (5, (10, 3)),  # long dash with offset
    (0, (5, 10)),  # loosely dashed
    'dashed',
    (0, (5, 1)),  # densely dashed
    (0, (3, 10, 1, 10)),  # loosely dash-dotted
    'dashdot',
    (0, (3, 1, 1, 1)),  # densely dash-dotted
    (0, (3, 5, 1, 5, 1, 5)),  # dash-dot-dotted
    (0, (3, 10, 1, 10, 1, 10)),  # loosely dash-dot-dotted
    (0, (3, 1, 1, 1, 1, 1)),  # densely dash-dot-dotted
]

# MY_MPL_STYLE = "$PYTHONPATH/feupy/visualization/styles/mystyle.mplstyle"

