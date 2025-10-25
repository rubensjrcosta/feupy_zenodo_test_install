"""
““” This script converts PSRCAT data to FITS format by querying the ATNF
pulsar database and processing specific pulsar parameters before saving
them in a FITS file. ““”

"""

from psrqpy import QueryATNF
from astropy.table import Table, Column
import logging

def format_table(table):
    """Format table"""
    for column in table.colnames:
        if column.startswith(("DEJ", "RAJ", "DIST")):
            table[column].format = ".3f"
        elif column.startswith(
            ("AGE", "EDOT", "BSURF", "P0", "P0_ERR")
        ):
            table[column].format = ".2e"
        elif column.startswith(
            ("sed_d")
        ):
            table[column].format = ".3e"

    return table

# Set up logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Define the pulsar parameters to query from ATNF
PSR_PARAMS = [
    'JNAME', 
    'RAJD', 
    'DECJD',
    'DIST',
    'DIST_DM', 
    'AGE', 
    'P0',
    'BSURF',
    'EDOT', 
    'TYPE', 
    'ASSOC',
]

try:
    # Query the ATNF pulsar catalog with the specified parameters
    log.info("Querying the ATNF pulsar catalog...")
    table = QueryATNF(params=PSR_PARAMS).table

    # Clean the TYPE column (replace '--' with None and ensure consistent string length)
    log.info("Processing TYPE column...")
    max_type_length = max(len(str(typ)) for typ in table['TYPE'] if typ != '--')  # Max length of valid 'TYPE'
    table['TYPE'] = [None if str(typ) == '--' else str(typ).ljust(max_type_length) for typ in table['TYPE']]
    table['TYPE'] = table['TYPE'].astype(f'<U{max_type_length}')  # Consistent dtype for the 'TYPE' column

    # Clean the ASSOC column (replace '--' with None)
    log.info("Processing ASSOC column...")
    table['ASSOC'] = [None if str(assoc) == '--' else str(assoc) for assoc in table['ASSOC']]
    max_assoc_length = max(len(str(assoc)) for assoc in table['ASSOC'] if assoc)
    table['ASSOC'] = table['ASSOC'].astype(f'<U{max_assoc_length}')  # Ensure consistent dtype for the 'ASSOC' column

    # Modify JNAME column to prepend 'PSR ' and ensure consistent string length
    log.info("Processing JNAME column...")
    table['NAME'] = [f"PSR {jname}" for jname in table['JNAME']]
    max_jname_length = max(len(str(jname)) for jname in table['NAME'])
    table['NAME'] = table['NAME'].astype(f'<U{max_jname_length}')  # Consistent dtype for 'JNAME'

    # Rename columns for RA and Dec with more descriptive names
    log.info("Renaming RAJD and DECJD columns...")
    table.rename_columns(['RAJD', 'DECJD'], ['RAJ2000', 'DEJ2000'])
    table.rename_columns(['RAJD_ERR', 'DECJD_ERR'], ['RAJ2000_ERR', 'DEJ2000_ERR'])

    # Update the metadata of the table
    data_name = "psrcat"
    table.meta['catalog_name'] = data_name

    # Write the table to a FITS file
    log.info(f"Writing data to {data_name}.fits...")
    table = table[
        'NAME',
        'JNAME',
        'DEJ2000',
        'DEJ2000_ERR',
        'RAJ2000',
        'RAJ2000_ERR',
        'DIST',
        'DIST_DM',
        'AGE',
        'EDOT',
        'BSURF',
        'P0',
        'P0_ERR',
        'ASSOC',
        'TYPE',
    ]
    table = format_table(table)  
    table.write(f"{data_name}_catalog.fits", format='fits', overwrite=True)

    log.info("FITS file successfully created.")

except Exception as e:
    log.error(f"An error occurred: {e}")

# from feupy.scripts.ipynb_to_gallery import convert_ipynb_to_gallery

# convert_ipynb_to_gallery('Untitled.ipynb', "make_psrcat.py")


# !pyflakes make_psrcat.py

