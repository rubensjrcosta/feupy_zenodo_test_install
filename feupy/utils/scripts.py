# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities to create scripts and command-line tools."""

from gammapy.utils.scripts import make_path
from astropy.coordinates import SkyCoord
import pickle


__all__ = [
    "is_documented_by",
    "pickling",
    "unpickling",
]


def is_documented_by(original):
    """
    Copy the docstring from the `original` function or class to a target.
    
    Parameters
    ----------
    original : function or class
        The function or class whose docstring is copied.

    Returns
    -------
    wrapper : function
        A decorator function that applies the original docstring to the target.
    """
    def wrapper(target):
        doc = '*** Docstring of internal function/class ***\n'
        if isinstance(original, list):
            for item in original:
                doc += f"{item.__qualname__}:\n{item.__doc__}\n"
        else:
            doc += f"{original.__doc__}\n"
        
        if target.__doc__:
            doc += f'\n*** Docstring of {target.__qualname__} ***\n{target.__doc__}'
        target.__doc__ = doc
        return target
    
    return wrapper

def pickling(object_instance, file_name):        
    """
    Serialize an object to a pickle file.
    
    Parameters
    ----------
    object_instance : object
        The object to serialize.
    file_name : str
        The name of the file (without extension) to store the pickle.
    """
    with open(make_path(f"{file_name}.pkl"), "wb") as fp:  
        pickle.dump(object_instance, fp)

def unpickling(file_name):        
    """
    Load an object from a pickle file.
    
    Parameters
    ----------
    file_name : str
        The name of the pickle file (without extension) to load.
    
    Returns
    -------
    object
        The deserialized object.
    """
    with open(make_path(f"{file_name}.pkl"), "rb") as fp:  
        return pickle.load(fp)
