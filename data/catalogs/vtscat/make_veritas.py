""" This script converts `<tev-source_id>.yaml` files to ECSV format
table.

VTSCat is the VERITAS Catalog of gamma-ray observations from this paper:
https://iopscience.iop.org/article/10.3847/2515-5172/acb147

"""

from feupy.scripts.ipynb_to_gallery import convert_ipynb_to_gallery


# Import necessary modules
import os
import logging
from astropy.table import Table, Column
import astropy.units as u

from feupy.utils.io import read_yaml
from gammapy.utils.table import table_row_to_dict
from gammapy.utils.scripts import make_path


# Set up logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Dictionary for source categories from VERITAS VTSCat
VTSCAT_DICT = {
    'detect': 'tev-0',
    'non_detect': 'tev-1',
    'crs': 'tev-3',
    'all': 'tev-'
}

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

def read_vtscat_sources(which='all', sources_path=SOURCES_PATH):
    """
    Read VERITAS sources from YAML files based on categories.

    Parameters
    ----------
    which : str, optional
        Key to select source type from `VTSCAT_DICT`, by default 'all'.
    sources_path : str or Path, optional
        Path to the directory containing YAML files, by default `SOURCES_PATH`.

    Returns
    -------
    list of dict
        List of parsed YAML data from selected source files.

    Raises
    ------
    ValueError
        If the source type provided does not exist in the catalog.
    FileNotFoundError
        If no matching files are found.
    """
    data_list = []
    if which not in VTSCAT_DICT:
        raise ValueError(f"Invalid category '{which}'. Valid options are {list(VTSCAT_DICT.keys())}")
    
    prefix = VTSCAT_DICT[which]
    for filename in os.listdir(sources_path):
        if filename.startswith(prefix) and filename.endswith('.yaml'):
            file_path = os.path.join(sources_path, filename)
            data = read_yaml(file_path)
            data_list.append(data)
            log.info(f"Loaded source: {data.get('source_name', filename)} from {filename}")

    if not data_list:
        raise FileNotFoundError(f"No matching files found for key '{which}' in directory '{sources_path}'.")
    
    return data_list

def flatten_data(data):
    """
    Flatten the nested dictionary structure from VTSCat sources.

    Parameters
    ----------
    data : dict
        Dictionary containing the source data.

    Returns
    -------
    dict
        Flattened data dictionary with default values for missing fields.
    """
    return {
        'veritas_id': data.get('source_id', None),
        'common_name': data.get('common_name', 'Unknown'),
        'veritas_name': data.get('veritas_name', {}).get('name', 'Unknown'),
        'veritas_components': data.get('veritas_name', {}).get('components', 'Unknown'),
        'other_names': ', '.join(data.get('other_names', ['Unknown'])),
        'where': data.get('where', 'Unknown'),
        'simbad_id': data.get('pos', {}).get('simbad_id', 'Unknown'),
        'ra': data.get('pos', {}).get('ra', None),
        'dec': data.get('pos', {}).get('dec', None),
        'reference_id': ', '.join(data.get('reference_id', ['Unknown']))
    }

# Read and process the data
data_list = read_vtscat_sources()
flattened_data_list = [flatten_data(data) for data in data_list]

# Define columns for the table based on the flattened data
columns = {
    'veritas_id': Column([d['veritas_id'] for d in flattened_data_list], dtype='int'),
    'common_name': Column([d['common_name'] for d in flattened_data_list], dtype='str', length=20),
    'veritas_name': Column([d['veritas_name'] for d in flattened_data_list], dtype='str', length=20),
    'veritas_components': Column([d['veritas_components'] for d in flattened_data_list], dtype='str', length=20),
    'other_names': Column([d['other_names'] for d in flattened_data_list], dtype='str', length=50),
    'where': Column([d['where'] for d in flattened_data_list], dtype='str', length=10),
    'simbad_id': Column([d['simbad_id'] for d in flattened_data_list], dtype='str', length=20),
    'ra': Column([d['ra'] for d in flattened_data_list], unit=u.deg, dtype='float'),
    'dec': Column([d['dec'] for d in flattened_data_list], unit=u.deg, dtype='float'),
    'reference_id': Column([d['reference_id'] for d in flattened_data_list], dtype='str', length=50)
}

# Create the Astropy Table and set metadata
table = Table(columns)
table.meta['catalog_name'] = 'vtscat'
table.meta['comments'] = ['Reference: https://zenodo.org/records/6163391']

# Set descriptions and formats for columns
table.sort('veritas_id')
table['ra'].description = 'Right Ascension (J2000)'
table['ra'].format = '.3f'
table['dec'].description = 'Declination (J2000)'
table['dec'].format = '.3f'
table['veritas_id'].description = 'Unique identifier of an object (integer)'
table['common_name'].description = 'Common name as frequently used in literature'
table['veritas_name'].description = 'VERITAS name as registered in the AU Registry'
table['other_names'].description = 'Other names used in literature or from catalogs (not a complete list)'
table['where'].description = 'Galactic (gal) or extragalactic (egal) object'
table['reference_id'].description = 'List of references for this catalog'
# Additional columns based on source data
table['source_name'] = [d['veritas_name'] if d['veritas_name'] != 'Unknown' else d['common_name'] for d in table]
table['source_name'].description = 'Source name: "veritas_name" if available, else "common_name"'
table['source_id'] = [_ for d in table]

# Populate the 'type' column based on veritas_id
table['type'] = ['detected' if d['veritas_id'] < 100000 else 
                         'non_detected' if d['veritas_id'] < 300000 else 'crs' 
                         for d in table]
table['type'].description = 'Type of the VERITAS type'
table = table['source_name', 'ra', 'dec', 'type', 'where', 'veritas_name', 'veritas_components', 'veritas_id','common_name', 'other_names', 'simbad_id', 'reference_id',]
# Save the table as ECSV format
output_path = f"{SOURCES_PATH}/vtscat.ecsv"
table.write(output_path, format='ascii.ecsv', overwrite=True)
log.info(f"Table saved as ECSV: {output_path}")

# convert_ipynb_to_gallery('Untitled.ipynb', "make_veritas.py")

# !pyflakes make_veritas.py

