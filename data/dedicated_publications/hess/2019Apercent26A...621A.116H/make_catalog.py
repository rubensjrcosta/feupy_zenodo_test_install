# Imports
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

def format_table(table):
    """Format table columns for better readability."""
    for column in table.colnames:
        if column.startswith(("Flux", 'sed_dnde')):
            table[column].format = ".3e"
        elif column.startswith(("ra", "dec", "e_min", "e_max", "e_ref", "sed_e", "sqrt_ts", "norm", "ts", "stat")):
            table[column].format = ".3f"
    return table

# Constants
BIBCODE = '2019A&A...621A.116H'
TAG = 'hess-2019A&A'
SED_TYPE = 'e2dnde'
REF_URL = "https://doi.org/10.1051/0004-6361/201834335"

DATA_FILES = {
    'HESS J1825-137': 'HESS_J1825-137.fits',
}

CATALOG = {
    "HESS J1825-137": {
        'position': SkyCoord('18 25 49 -13 46 35', unit=(u.hourangle, u.deg)),  # HESS Collaboration (2019)
        'radius': 0.8*u.deg,
        'radius_core': 0.4*u.deg,
        'spec_type': 'ecpl',
    },
}

# Prepare Data for Table
source_names = []
ra_list = []
dec_list = []
radius = []
radius_core = []
spec_model_type = []

catalog_table_meta = {
    'catalog_name': TAG,
    "SED_TYPE": SED_TYPE,
    "reference": REF_URL,
}

# Generate Table Columns
for source_name, source_data in CATALOG.items():
    source_names.append(source_name)
    ra_list.append(source_data['position'].ra.deg)
    dec_list.append(source_data['position'].dec.deg)
    radius.append(source_data['radius'])
    radius_core.append(source_data['radius_core'])
    spec_model_type.append(source_data['spec_type'])

# Create and Save the Table
table = Table(
    [source_names, ra_list, dec_list, radius, radius_core, spec_model_type],
    names=('source_name', 'ra', 'dec', 'radius', 'radius_core', 'spec_type')
)
table.meta.update(catalog_table_meta)
table['ra'].unit = 'deg'
table['dec'].unit = 'deg'
table = format_table(table)
table.write('catalog.ecsv', format='ascii.ecsv', overwrite=True)
