# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Sources class."""

import collections.abc
import copy
import yaml
import numpy as np
import logging
from astropy.coordinates import SkyCoord
from feupy.catalog import FEUPY_CATALOG_REGISTRY

log = logging.getLogger(__name__)


__all__ = [
    "Sources",
    "get_catalog_tag",
]


def get_catalog_tag(source):
    """
    Retrieve the catalog tag for a given source.

    This function checks if the input source is an instance of any class 
    in the `CATALOG_REGISTRY` and returns the corresponding catalog tag.

    Parameters
    ----------
    source : object
        The source object for which the catalog tag is to be retrieved.
        The source should be an instance of one of the classes in 
        `CATALOG_REGISTRY`.

    Returns
    -------
    str
        The tag of the first matching catalog from `FEUPY_CATALOG_REGISTRY`.

    Raises
    ------
    ValueError
        If the input `source` does not match any class in `FEUPY_CATALOG_REGISTRY`, 
        a `ValueError` is raised indicating no matching catalog was found.

    Examples
    --------
    Example usage of the `get_catalog_tag` function:

    >>> source = SourceCatalogObjectGammaCat()
    >>> tag = get_catalog_tag(source)
    >>> print(tag)
    'gamma-cat'

    Notes
    -----
    This function assumes that `CATALOG_REGISTRY` is a list of catalog 
    classes and that each class in the registry has a `tag` attribute. 
    The function returns the tag of the first catalog that matches 
    the type of the given source.
    """
    # Find the first matching catalog
    matching_catalog = next(
        (catalog for catalog in FEUPY_CATALOG_REGISTRY if isinstance(source, catalog.source_object_class)), 
        None
    )

    if matching_catalog is None:
        log.error(f"Failed to found catalog for source: {source}")
        raise ValueError(f"No matching catalog found for source: {source}")

    return matching_catalog.tag


class Sources(collections.abc.MutableSequence):
    """Source collection.

    Parameters
    ----------
    datasets : `Source` or list of `Source`
        Sources
    """

    def __init__(self, sources=None):
        if sources is None:
            sources = []
        elif isinstance(sources, Sources):
            sources = sources._sources
        elif any(isinstance(sources, _.source_object_class) for _ in FEUPY_CATALOG_REGISTRY):
            sources = [sources]
        elif not isinstance(sources, list):
            log.error(f"Failed Invalid type: {sources!r}")
            raise TypeError(f"Invalid type: {sources!r}")

        unique_labels = []
        for source in sources:
            label = self._set_source_label(source)
            if label in unique_labels:
                log.error(f"Failed Source name '{source.name}' from {get_catalog_tag(source)} already exists!")
                raise ValueError(f"Source name '{source.name}' from {get_catalog_tag(source)} already exists!")
            unique_labels.append(label)

        self._sources = sources

    def __getitem__(self, key):
        return self._sources[self.index(key)]

    def __delitem__(self, key):
        del self._sources[self.index(key)]

    def __setitem__(self, key, source):
        if any(isinstance(source, _) for _ in self._sources):
            label = self._set_source_label(source)
#             label = (source_name, catalog_tag)
            if label in self.labels:
                log.error(f"Failed Source name '{source.name}' from {get_catalog_tag(source)} already exists!")
                raise ValueError(f"Source name '{source.name}' from {get_catalog_tag(source)} already exists!")
            self._sources[self.index(key)] = source
        else:
            log.error(f"Failed Invalid type: {type(source)!r}")
            raise TypeError(f"Invalid type: {type(source)!r}")

    def insert(self, idx, source):
        if any(isinstance(source, _.source_object_class) for _ in FEUPY_CATALOG_REGISTRY):
            label = self._set_source_label(source)
            if label in self.labels:
                log.error(f"Failed Source name '{source.name}' from {get_catalog_tag(source)} already exists!")
                raise ValueError(f"Source name '{source.name}' from {get_catalog_tag(source)} already exists!")
            self._sources.insert(idx, source)
        else:
            log.error(f"Failed Invalid type: {type(source)!r}")
            raise TypeError(f"Invalid type: {type(source)!r}")

    def _set_source_label(cls, source):
        tag = get_catalog_tag(source)
        return f'{source.name} ({tag})'
    
    def index(self, key):
        if isinstance(key, (int, slice)):
            return key
        elif isinstance(key, str):
            return self.names.index(key)
        elif any(isinstance(key, _) for _ in self._sources):
            return self._sources.index(key)
        else:
            log.error(f"Failed Invalid type: {type(key)!r}")
            raise TypeError(f"Invalid type: {type(key)!r}")
       
    def __len__(self):
        return len(self._sources)
    
    def copy(self):
        """A deep copy."""
        return copy.deepcopy(self)
    
    @property
    def names(self):
        return [_.name for _ in self._sources]

    @property
    def labels(self):
        return [self._set_source_label(_) for _ in self._sources]
    
    @property
    def positions(self):
        """Source positions as a `~astropy.coordinates.SkyCoord` object."""
        ra = [_.position.icrs.ra for _ in self._sources]
        dec = [_.position.icrs.dec for _ in self._sources]
        return SkyCoord(ra, dec, frame='icrs')
    
    def write(self, filename, overwrite=False):
        _dict = {}
        dict_sources = {}

        for index, source in enumerate(self._sources):
            _dict[index] = dict(
                name=source.name,
                catalog=get_catalog_tag(source),
            )

        dict_sources['Sources'] = _dict

        yaml_string = yaml.dump(dict_sources)
        sources = yaml.safe_load(yaml_string)

        with open(filename, 'w') as file:
            yaml.dump(sources, file)
            
    def read(self, filename):
        with open(filename, 'r') as stream:
            try:
                d = yaml.safe_load(stream)
                self._sources = [
                    FEUPY_CATALOG_REGISTRY.get_cls(d['Sources'][_]['catalog'])()[d['Sources'][_]['name']]  
                    for _ in d['Sources'].keys()
                ]
            except yaml.YAMLError as error:
                print(error)