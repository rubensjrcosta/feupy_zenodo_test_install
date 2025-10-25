"""
“““This script converts extraHAWC data to FITS format.”“”

"""

from astropy.table import Table, Column, MaskedColumn
from astropy.coordinates import SkyCoord
import astropy.units as u
from feupy.utils.table import pad_list_to_length
from gammapy.utils.scripts import make_path

import numpy as np
from gammapy.utils.scripts import make_path
from gammapy.modeling.models import LogParabolaSpectralModel, ExpCutoffPowerLawSpectralModel, PowerLawSpectralModel

from feupy.utils.units import FRAME_ICRS, UNIT_DEG, DEFAULT_ENERGY_UNIT, DEFAULT_SED_UNIT


def read_sed_data(file, length):
    """Read spectral energy distribution (SED) data from a file."""
    e_ref_list, e2dnde_list, e2dnde_errn_list, e2dnde_errp_list, e2dnde_ul_list, is_ul_list = [], [], [], [], [], []
    table = Table.read(file) 
    return (
        table['e_ref'],
        table['e2dnde'],
        table['e2dnde_errn'],
        table['e2dnde_errp'],
        table['e2dnde_ul'],
        table['is_ul'],
    )

def format_table(table):
    """Format table columns for better readability."""
    for column in table.colnames:
        if column.startswith(('sed_e2dnde')):
            table[column].format = ".3e"
        elif column.startswith(("e_min", "e_max", "e_ref", "sed_e", "sqrt_ts", "norm", "ts", "stat")):
            table[column].format = ".3f"
    return table

# Constants
BIBCODE = '2021ApJ...907L..30A'
TAG = 'hawc-2021ApJ'
SED_TYPE = 'e2dnde'
DEFAULT_LENGTH = 10
REFERENCE_ENERGY = 18 * u.TeV
REF_URL = "https://iopscience.iop.org/article/10.3847/2041-8213/abd77b"


DATA_FILES = {
    'HAWC J1825-138': 'HAWC_J1825-138.fits',
    'HAWC J1826-128': 'HAWC_J1826-128.fits',
    'HAWC J1825-134': 'HAWC_J1825-134.fits',
}

CATALOG = {
    "HAWC J1825-138": {
        "position":  SkyCoord(276.38, -13.86, unit='deg', frame='icrs'),
        'spectral_model': ExpCutoffPowerLawSpectralModel(
            amplitude=2.7e-14 * u.Unit("cm-2 s-1 TeV-1"),
            index=2.02,
            lambda_=1/27 * u.Unit("TeV-1"),
            reference=18 * u.TeV,),
        },
    "HAWC J1826-128": {
        "position":  SkyCoord(276.50, -12.86, unit='deg', frame='icrs'),
        'spectral_model': ExpCutoffPowerLawSpectralModel(
            amplitude=2.7e-14 * u.Unit("cm-2 s-1 TeV-1"),
            index=1.2,
            lambda_=(1/24) * u.Unit("TeV-1"),
            reference=18 * u.TeV,),
        },
    "HAWC J1825-134": {
        "position":  SkyCoord(276.44, -13.42, unit='deg', frame='icrs'),
        'spectral_model': PowerLawSpectralModel(
            index=2.28,
            amplitude="4.2e-15 TeV-1 cm-2 s-1",
            reference=18 * u.TeV,),
        },
}

# Data Preparation

# Prepare data for the table
source_names = []
ra_list = []
dec_list = []
spec_model_type = []

catalog_table_meta = {'catalog_name': TAG,
                      "SED_TYPE": SED_TYPE,
                      "reference": REF_URL
}

# Generate Table Columns
# Iterate over the sources to extract values
for source_name, source_data in CATALOG.items():
    source_names.append(source_name)
    
    ra_list.append(source_data['position'].ra.deg)
    dec_list.append(source_data['position'].dec.deg)
    
    # Extract spectral model details
    model = source_data['spectral_model']
    spec_model_type.append(model.tag[1])
        
table = Table([source_names, ra_list, dec_list,spec_model_type],
              names=('source_name', 'ra', 'dec', 'spec_type'))
table.meta.update(catalog_table_meta)

table.write('catalog.ecsv', format='ascii.ecsv', overwrite=True)