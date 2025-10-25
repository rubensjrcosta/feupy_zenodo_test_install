# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Script to collect pulsar spectral flux data from a catalogue,
convert NumPy objects to plain Python types, and save as a YAML file.
"""

import yaml
import numpy as np
from pulsar_spectra.catalogue import collect_catalogue_fluxes


def convert_numpy_objects(obj):
    """
    Recursively convert NumPy objects to plain Python types.
    
    Parameters
    ----------
    obj : object
        The object to be converted. Can be a dictionary, list, tuple, NumPy array, or scalar.
    
    Returns
    -------
    object
        The converted object with all NumPy-specific types replaced by plain Python types.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.generic):
        return obj.item()  # Convert NumPy scalar to Python scalar
    elif isinstance(obj, dict):
        return {k: convert_numpy_objects(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_objects(x) for x in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_numpy_objects(x) for x in obj)
    return obj


def main(output_file="pulsar_spectra.yaml"):
    """
    Main function to collect pulsar catalogue data and save it as a YAML file.
    
    Parameters
    ----------
    output_file : str, optional
        Name of the output YAML file. Default is 'pulsar_spectra.yaml'.
    """
    # Collect the pulsar catalogue data
    cat_dict = collect_catalogue_fluxes()
    
    # Convert NumPy objects to plain Python types for YAML compatibility
    safe_cat_dict = convert_numpy_objects(cat_dict)
    
    # Save the dictionary to a YAML file
    with open(output_file, "w") as yaml_file:
        yaml.dump(safe_cat_dict, yaml_file, default_flow_style=False)
    print(f"Pulsar spectral data saved to {output_file}")


if __name__ == "__main__":
    main()

