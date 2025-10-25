# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for table manipulation.

This module provides functions for handling and processing tables, including:
    - Padding lists for table column compatibility.
    - Removing rows with NaN values.
    - Reading and writing tables in CSV and FITS formats.
"""

import sys, os
import numpy as np
from astropy.table import Table
from gammapy.utils.scripts import make_path

__all__ = [
    "pad_list_to_length",
    "remove_nan_rows",
    "write_table",
    "read_table",
]





def pad_list_to_length(length, input_list):
    """
    Pads a list with `None` values until it reaches the specified length.
    If the list is already longer than the target length, an AttributeError is raised.

    Parameters
    ----------
    length : int
        The desired length of the list after padding.
    input_list : list
        The list to be padded with `None` values.

    Returns
    -------
    list
        A list with length equal to `length`, padded with `None` values if necessary.

    Raises
    ------
    AttributeError
        If the `input_list` is already longer than `length`.
    """
    diff_len = length - len(input_list)
    if diff_len < 0:
        raise AttributeError('The input list is longer than the specified target length.')
    return input_list + [np.nan] * diff_len

def remove_nan_rows(table):
    """
    Remove rows containing NaN values from an Astropy table.

    Parameters
    ----------
    table : astropy.table.Table
        The input Astropy table from which rows containing NaN values will be removed.

    Returns
    -------
    astropy.table.Table
        A new table with all rows containing NaN values removed.
    """
    has_nan = np.zeros(len(table), dtype=bool)
    for col in table.itercols():
        if col.info.dtype.kind == 'f':  # Check if the column is of floating-point type
            has_nan |= np.isnan(col)
    return table[~has_nan]

def write_table(table, path_file, file_name, overwrite=False):
    """
    Write an Astropy table to a file in CSV or FITS format, determined by the file extension.
    
    Parameters
    ----------
    table : astropy.table.Table
        Table to write to a file.
    path_file : str
        Path to the directory where the file will be saved.
    file_name : str
        Name of the file, with `.csv` or `.fits` extension to specify format.
    overwrite : bool, optional
        Whether to overwrite the file if it already exists (default is False).

    Raises
    ------
    ValueError
        If the file extension is not `.csv` or `.fits`.
    """
    # Determine file format based on extension
    if file_name.endswith(".csv"):
        file_format = "ascii.ecsv"
    elif file_name.endswith(".fits"):
        file_format = "fits"
    else:
        raise ValueError("File name must end with .csv or .fits")

    # Construct the full path
    path_os = os.path.abspath(os.path.join(path_file, file_name))

    # Write the table to the specified format
    table.write(path_os, format=file_format, overwrite=overwrite)


def read_table(path_file, file_name):
    """
    Read an Astropy table from a file in CSV or FITS format, determined by the file extension.
    
    Parameters
    ----------
    path_file : str
        Path to the directory where the file is located.
    file_name : str
        Name of the file, with `.csv` or `.fits` extension to specify format.
    
    Returns
    -------
    astropy.table.Table
        The table read from the specified file.
    
    Raises
    ------
    ValueError
        If the file extension is not `.csv` or `.fits`.
    """
    # Determine file format based on extension
    if file_name.endswith(".csv"):
        file_format = "ascii.ecsv"
    elif file_name.endswith(".fits"):
        file_format = "fits"
    else:
        raise ValueError("File name must end with .csv or .fits")

    # Construct the full path
    path_os = os.path.abspath(os.path.join(path_file, file_name))

    # Read the table from the specified format
    return Table.read(path_os, format=file_format)

def write_tables_fits(table, path_file, file_name, overwrite = True):
    # Writes the flux points table in the fits format
    path_os = os.path.abspath(
        os.path.join(
            f"{path_file}/{file_name}.fits"
        )
    )      

    if path_os not in sys.path:
        sys.path.append(path_os)

    table.write(
        f"{path_os}",
        format = 'fits', 
        overwrite = overwrite
    )   
    return

def write_tables_csv(table, path_file, file_name, overwrite = True):
# Writes the flux points table in the csv format
    path_os = os.path.abspath(
        os.path.join(
            f"{path_file}/{file_name}.csv"
        )
    )

    if path_os not in sys.path:
        sys.path.append(path_os)

    table.write(
        f"{path_os}",
        format = 'ascii.ecsv', 
        overwrite = overwrite
    )   
    return

# In[1]:


# from astropy.units import Quantity

# def table_row_to_dict(row, make_quantity=True):
#     """Make one source data dictionary.

#     Parameters
#     ----------
#     row : `~astropy.table.Row`
#         Row.
#     make_quantity : bool, optional
#         Make quantity values for columns with units.
#         Default is True.

#     Returns
#     -------
#     data : dict
#         Row data.
#     """
    
#     data = {}
#     for name, col in row.columns.items():
#         val = row[name]

#         if make_quantity and col.unit:
#             val = Quantity(val, unit=col.unit)
#         data[name] = val
#     return data


# In[ ]:


# # Data column for use in a Table object to string using join()
# def column_to_string(column):
#     return f'[{",".join(str(element) for element in list(column))}]'


# # In[ ]:


# def append_nones(length, list_):
#     """
#     Appends Nones to list to get length of list equal to `length`.
#     If list is too long raise AttributeError
#     """
#     diff_len = length - len(list_)
#     if diff_len < 0:
#         raise AttributeError('Length error list is too long.')
#     return list_ + [None] * diff_len


# # In[ ]:


# def remove_nan(table):
#     """Remove lines containing a nan"""
#     has_nan = np.zeros(len(table), dtype=bool)
#     for col in table.itercols():
#         if col.info.dtype.kind == 'f':
#             has_nan |= np.isnan(col)
#     return table[~has_nan]

