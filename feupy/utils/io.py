# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities for working with I/O."""

import yaml
import logging
from pathlib import Path

# Set up logging as per Gammapy's development guidelines
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def mkdir_sub_directory(parent_directory, child_directory=None):
    """Create parent and optional child directories, and return their paths.

    Parameters
    ----------
    parent_directory : str or Path
        Path to the parent directory.
    child_directory : str, optional
        Name of the child directory, by default None.

    Returns
    -------
    Path or tuple of Paths
        If only `parent_directory` is provided, returns Path to the parent directory.
        If both `parent_directory` and `child_directory` are provided, returns a tuple with Paths to parent and child directories.

    Examples
    --------
    >>> parent_dir = mkdir_sub_directory('data')
    >>> parent_dir, child_dir = mkdir_sub_directory('data', 'subfolder')
    """
    path_parent = Path(parent_directory)
    path_parent.mkdir(parents=True, exist_ok=True)
    log.info(f"Directory '{path_parent}' created.")

    if child_directory is not None:
        path_child = path_parent / child_directory
        path_child.mkdir(parents=True, exist_ok=True)
        log.info(f"Directory '{path_child}' created.")
        return path_parent, path_child

    return path_parent


def read_yaml(file_path):
    """Read and parse a YAML file.

    Parameters
    ----------
    file_path : str or Path
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed YAML content as a dictionary.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    yaml.YAMLError
        If there is an error parsing the YAML file.

    Examples
    --------
    >>> data = read_yaml('data/info.yaml')
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        log.error(f"File {file_path} not found.")
        raise FileNotFoundError(f"File {file_path} not found.")

    try:
        with file_path.open('r') as file:
            data = yaml.safe_load(file)
            log.info(f"Successfully loaded YAML file: {file_path}")
            return data
    except yaml.YAMLError as exc:
        log.error(f"Error parsing YAML file {file_path}: {exc}")
        raise
