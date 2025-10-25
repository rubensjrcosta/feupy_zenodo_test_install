# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""String handling utility functions."""

import json

def string_to_filename_format(name):
    """
    Convert a name string to a filename-friendly format.

    This function replaces spaces with underscores and removes colons from
    the input string, making it suitable for use as a filename.

    Parameters
    ----------
    name : str
        The input string, typically a source or object name.

    Returns
    -------
    str
        The modified string suitable for use as a filename.

    Examples
    --------
    >>> string_to_filename_format("HESS J1825:137")
    'HESS_J1825137'

    >>> string_to_filename_format("Crab Nebula")
    'Crab_Nebula'
    """
    return name.replace(" ", "_").replace(":", "")


def string_to_list(string):
    """
    Convert a string representation of a list into an actual Python list.

    This function uses `json.loads()` to interpret a string containing
    a JSON-like list and returns it as a Python list.

    Parameters
    ----------
    string : str
        A string representing a list in JSON format.

    Returns
    -------
    list
        A Python list derived from the string input.

    Examples
    --------
    >>> string_to_list("[1, 2, 3]")
    [1, 2, 3]

    >>> string_to_list('["a", "b", "c"]')
    ['a', 'b', 'c']
    """
    return json.loads(string)






