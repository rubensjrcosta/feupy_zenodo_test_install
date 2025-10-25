# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Configuration and utilities for 1LHAASO analysis."""

# Imports
from gammapy.utils.scripts import make_path

# ------------------------------------------------------------------------------
# Paths and directories for data, datasets, and figures
# ------------------------------------------------------------------------------

DATA_PATH = '/home/born-again/Documents/GitHub/feupy-dev/data/vtscat/'
DATA_PATH = make_path(DATA_PATH)

DATA_PATH.mkdir(parents=True, exist_ok=True)

DATASETS_PATH = make_path(f'{DATA_PATH}/datasets/')
DATASETS_PATH.mkdir(parents=True, exist_ok=True)

SOURCES_PATH = make_path(f'{DATA_PATH}/sources/')
SOURCES_PATH.mkdir(parents=True, exist_ok=True)
