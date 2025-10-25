# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides functions to read and write reference markers in YAML format.

Functions:
    - write_ref_markers: Writes a dictionary of reference markers to a YAML file.
    - read_ref_markers: Reads reference markers from a YAML file and returns them as a dictionary.
"""

import random
import yaml
import os 
from gammapy.utils.scripts import recursive_merge_dicts
from gammapy.datasets import Datasets
from feupy.catalog import CATALOG_REGISTRY
from feupy.visualization import PALETTE_DEFAULT, MARKERS_DEFAULT, MARKERS_DEFAULT_DICT, LINESTYLES_DEFAULT
from feupy.sources import Sources
from feupy.sources import get_catalog_tag
# from feupy.scripts.ipynb_to_gallery import convert_ipynb_to_gallery
import random
import astropy.units as u



__all__ = ["write_ref_markers", "read_ref_markers", "generate_specified_marker_set", "generate_catalog_markers", "get_colors", "get_linestyles", "get_markers", "get_kwargs_fit" ]


def generate_repeated_list(items, num_items, shuffle=False):
    """
    Generate a list with repeated elements from a given list of items.

    Parameters
    ----------
    items : list
        List of items (e.g., colors, markers, linestyles).
    num_items : int
        Number of items to generate.
    shuffle : bool, optional
        If True, shuffle the items list each time it's repeated.

    Returns
    -------
    list
        A list containing the requested number of items, repeated if needed.
    """
    result = (items * (num_items // len(items))) + items[:num_items % len(items)]
    if shuffle:
        random.shuffle(result)
    return result


def get_colors(palette, num_colors, shuffle=False):
    """
    Generate a list of colors from the given palette.

    Parameters
    ----------
    palette : list
        List of available colors.
    num_colors : int
        Number of colors to generate.
    shuffle : bool, optional
        If True, shuffle the colors list.

    Returns
    -------
    list
        A list of colors.
    """
    return generate_repeated_list(palette, num_colors, shuffle)


def get_markers(markers, num_markers, shuffle=False):
    """
    Generate a list of markers from the given marker set.

    Parameters
    ----------
    markers : list
        List of available markers.
    num_markers : int
        Number of markers to generate.
    shuffle : bool, optional
        If True, shuffle the markers list.

    Returns
    -------
    list
        A list of markers.
    """
    return generate_repeated_list(markers, num_markers, shuffle)


def get_linestyles(linestyles, num_linestyles, shuffle=False):
    """
    Generate a list of linestyles from the given styles.

    Parameters
    ----------
    linestyles : list
        List of available linestyles.
    num_linestyles : int
        Number of linestyles to generate.
    shuffle : bool, optional
        If True, shuffle the linestyles list.

    Returns
    -------
    list
        A list of linestyles.
    """
    return generate_repeated_list(linestyles, num_linestyles, shuffle)


def write_ref_markers(ref_markers, filename, overwrite=False):
    """
    Write a dictionary of reference markers to a YAML file.
    
    Parameters
    ----------
    ref_markers : dict
        Dictionary containing reference markers to write to the file.
    filename : str
        Name of the YAML file where data will be written.
    overwrite : bool, optional
        Whether to overwrite the file if it already exists (default is False).
        
    Raises
    ------
    FileExistsError
        If `overwrite` is False and the file already exists.
    """
    # Prepare dictionary for YAML dumping
    dict_ref_markers = {'ref_markers': ref_markers}

    # Check if file exists and overwrite is disabled
    if not overwrite and os.path.isfile(filename):
        raise FileExistsError(f"The file '{filename}' already exists. Set overwrite=True to overwrite.")

    # Write YAML data to file
    with open(filename, 'w') as file:
        yaml.dump(dict_ref_markers, file)


def read_ref_markers(filename):
    """
    Read reference markers from a YAML file and return them as a dictionary.
    
    Parameters
    ----------
    filename : str
        Name of the YAML file to read.
        
    Returns
    -------
    dict
        Dictionary of reference markers read from the file.
    
    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    yaml.YAMLError
        If there is an error while parsing the YAML file.
    """
    # Verify if file exists before attempting to read
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' does not exist.")

    # Read YAML data from file
    try:
        with open(filename, 'r') as stream:
            data = yaml.safe_load(stream)
            return data.get('ref_markers', {})
    except yaml.YAMLError as error:
        print(f"Error reading YAML file '{filename}': {error}")
        raise

        
def generate_specified_marker_set(refs, marker='o', marker_size=6, uniform_size=True, PALETTE=None):
    """
    Generate a set of markers with the specified marker and unique colors from PALETTE_DEFAULT.

    Parameters
    ----------
    refs : list
        List of references (e.g., dataset names).
    marker : str, optional
        Marker shape (default is 'o' for circle).
    marker_size : int, optional
        Size of the marker (default is 6).
    uniform_size : bool, optional
        If True, all markers will have uniform size, by default True.
    PALETTE : list or None, optional
        List of colors to choose from (default is PALETTE_DEFAULT).

    Returns
    -------
    dict
        Dictionary containing markers for each reference with unique colors.

    Raises
    ------
    ValueError
        If there are not enough unique colors in the palette to assign one to each reference.
    """
    if PALETTE is None:
        PALETTE = PALETTE_DEFAULT

    num_refs = len(refs)

    if num_refs > len(PALETTE):
        raise ValueError("Not enough unique colors in PALETTE to assign to each reference.")
    
    if uniform_size:
        marker_size = round(marker_size * MARKERS_DEFAULT_DICT[marker][0], 2)

    random.shuffle(PALETTE)

    ref_markers = {}
    for i, ref in enumerate(refs):
        color = PALETTE[i][0]
        ref_markers[ref] = {
            'label': ref,
            'marker': marker,
            'color': color,
            'markersize': marker_size,
        }

    return ref_markers

# # Define a function to map catalog tags to default markers
# def map_catalog_tags_to_markers(catalog_registry, markers_default):
#     """
#     Maps catalog tags from `catalog_registry` to marker symbols from `markers_default`.

#     Parameters:
#     - catalog_registry: A registry containing catalog objects, each with a 'tag' attribute.
#     - markers_default: A list of lists, where each inner list contains marker info as
#       [symbol, scale, description].

#     Returns:
#     - dict: A dictionary mapping catalog tags to marker symbols.
#     """
#     dict_cat = {}

#     # Extract tags from catalog registry
#     tags = [cat.tag for cat in catalog_registry]
    
#     # Map each tag to a marker symbol
#     for index, tag in enumerate(tags):
#         # Ensure that index is within bounds of markers_default
#         if index < len(markers_default):
#             dict_cat[tag] = markers_default[index][0]
#         else:
#             dict_cat[tag] = None  # Handle cases where there are more tags than markers

#     return dict_cat

# Define a function to map catalog tags to default markers
def map_catalog_tags_to_markers(sources):
    """
    Maps catalog tags from `catalog_registry` to marker symbols from `markers_default`.

    Parameters:
    - catalog_registry: A registry containing catalog objects, each with a 'tag' attribute.
    - markers_default: A list of lists, where each inner list contains marker info as
      [symbol, scale, description].

    Returns:
    - dict: A dictionary mapping catalog tags to marker symbols.
    """
    dict_cat = {}

    # Extract tags from catalog registry
    tags = [get_catalog_tag(source) for source in sources]
    for tag in tags:
        if 'fhl' in tag or 'fgl' in tag: 
            dict_cat[tag] = 'v'
        elif 'hgps' in tag or 'hess' in tag: 
            dict_cat[tag] = 'p'
        elif 'vtscat' in tag or 'veritas' in tag: 
            dict_cat[tag] = '>'
        elif 'hwc' in tag:
            dict_cat[tag] = '8'
        elif 'gamma' in tag:
            dict_cat[tag] = 'h'        
        elif 'LHAASO' in tag:
            dict_cat[tag] = 's'
        elif any(_ == tag for _ in ['psrcat', '2PC', '3PC']):
            dict_cat[tag] = '*'
        else:
            raise ValueError(f"No marker to catalog: {tag}.")

    return dict_cat



# def generate_catalog_markers(sources, datasets=None, marker_size=6, MARKERS=MARKERS_DEFAULT, CATALOG_REGISTRY=CATALOG_REGISTRY, PALETTE=None):
#     """
#     Generate a dictionary of markers for a given set of sources based on their catalog tags.

#     Parameters
#     ----------
#     sources : list of `Source` or Datasets
#         List of source objects. Each source should have a catalog tag retrievable via `get_catalog_tag`.
#     marker_size : int, optional
#         Size of the marker to be used for plotting, by default 6.
#     PALETTE : list or None, optional
#         Color palette for the markers, by default None.

#     Returns
#     -------
#     dict
#         Dictionary where the keys are source names and the values are dictionaries containing the marker configurations.

#     Raises
#     ------
#     ValueError
#         If `get_catalog_tag` does not match a known catalog for any source.
#     """
#     ref_markers = {}
    
#     catalog_markers = map_catalog_tags_to_markers(CATALOG_REGISTRY, MARKERS_DEFAULT)
    
#     for tag in catalog_markers.keys():
#         sources_selected = [source for source in sources if get_catalog_tag(source) == tag]
#         if sources_selected:
#             sources_selected = Sources(sources_selected)
#             labels = sources_selected.labels
#             if datasets:
#                 for source_selected in sources_selected:
#                     datasets_names = [_.name for _ in datasets if _.name.startswith(source_selected.name) ]
#                 if datasets_names:
#                     labels.extend(datasets_names)

#             _ref_markers = generate_specified_marker_set(
#                labels,
#                 marker=catalog_markers[tag],
#                 marker_size=marker_size,
#                 PALETTE=PALETTE
#             )
#             ref_markers = recursive_merge_dicts(ref_markers, _ref_markers)
#     return ref_markers
def generate_catalog_markers(sources, datasets=None, marker_size=6, MARKERS=MARKERS_DEFAULT, CATALOG_REGISTRY=CATALOG_REGISTRY, PALETTE=None):
    """
    Generate a dictionary of markers for a given set of sources based on their catalog tags.

    Parameters
    ----------
    sources : list of `Source` or Datasets
        List of source objects. Each source should have a catalog tag retrievable via `get_catalog_tag`.
    marker_size : int, optional
        Size of the marker to be used for plotting, by default 6.
    PALETTE : list or None, optional
        Color palette for the markers, by default None.

    Returns
    -------
    dict
        Dictionary where the keys are source names and the values are dictionaries containing the marker configurations.

    Raises
    ------
    ValueError
        If `get_catalog_tag` does not match a known catalog for any source.
    """
    ref_markers = {}
    
    catalog_markers = map_catalog_tags_to_markers(sources)
    
    for tag in catalog_markers.keys():
        marker = catalog_markers[tag]
        sources_selected = Sources([source for source in sources if catalog_markers[get_catalog_tag(source)] == marker])            
        labels = sources_selected.labels
        if datasets:
            for source_selected in sources_selected:
                datasets_names = [_.name for _ in datasets if _.name.startswith(source_selected.name) ]
            if datasets_names:
                labels.extend(datasets_names)

        _ref_markers = generate_specified_marker_set(
           labels,
            marker=catalog_markers[tag],
            marker_size=marker_size,
            PALETTE=PALETTE
        )
        ref_markers = recursive_merge_dicts(ref_markers, _ref_markers)
    return ref_markers

def get_kwargs_fit(refs, energy_bounds = [[5e-2, 2e3] * u.TeV], get_color=False, marker = ','):
    if len(energy_bounds) == len(refs):
        _energy_bounds = energy_bounds
    else:
        if len(energy_bounds) == 1:
            _energy_bounds = [energy_bounds[0]]*len(refs)
        else:
            raise ValueError("The length of the energy bounds list must be one or equal to the length of the references")
        
    kwargs_legend = {}
    linestyles = LINESTYLES_DEFAULT
    if len(linestyles) < len(refs):
        while len(linestyles) < len(refs)+1:
            linestyles.extend(linestyles)

    kwargs_legend = {}
    _color = 'black'
    for index, (label, ls, color) in enumerate(zip(refs, linestyles, PALETTE_DEFAULT)):
        if get_color:
            _color = color
        kwargs_legend[index] = dict(
            energy_bounds=_energy_bounds[index],
            label = label,
            ls = ls,
            marker=marker,
            color = _color,
            )
    return kwargs_legend
