# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Catalog utilities classes."""

import logging
import numpy as np
from typing import List, Optional
from feupy.catalog import FEUPY_CATALOG_REGISTRY

log = logging.getLogger(__name__)

catalog_2fhl = FEUPY_CATALOG_REGISTRY.get_cls('2fhl')()
catalog_3fhl = FEUPY_CATALOG_REGISTRY.get_cls('3fhl')()

catalog_3fgl = FEUPY_CATALOG_REGISTRY.get_cls('3fgl')()
catalog_4fgl = FEUPY_CATALOG_REGISTRY.get_cls('4fgl')()

catalog_2hwc = FEUPY_CATALOG_REGISTRY.get_cls('2hwc')()
catalog_3hwc = FEUPY_CATALOG_REGISTRY.get_cls('3hwc')()
catalog_ehwc = FEUPY_CATALOG_REGISTRY.get_cls('ehwc')()
catalog_extra_hawc = FEUPY_CATALOG_REGISTRY.get_cls('hwc-2021ApJ')()

catalog_hgps = FEUPY_CATALOG_REGISTRY.get_cls('hgps')()
catalog_extra_hess = FEUPY_CATALOG_REGISTRY.get_cls('hess-2019A&A')()

catalog_gamma_cat = FEUPY_CATALOG_REGISTRY.get_cls('gamma-cat')()

catalog_vtscat = FEUPY_CATALOG_REGISTRY.get_cls('vtscat')()

catalog_veritas = FEUPY_CATALOG_REGISTRY.get_cls('veritas-2018ApJ')()

catalog_lhaaso = FEUPY_CATALOG_REGISTRY.get_cls('LHAASO')()
catalog_1lhaaso = FEUPY_CATALOG_REGISTRY.get_cls('1LHAASO')()
catalog_extra_lhaaso = FEUPY_CATALOG_REGISTRY.get_cls('LHAASO-2024icrc')()

catalog_psrcat = FEUPY_CATALOG_REGISTRY.get_cls('psrcat')()

def load_catalogs(catalogs: Optional[List] = FEUPY_CATALOG_REGISTRY) -> List:
    """Load a list of catalogs from the provided registry.

    Parameters:
    -----------
    catalogs : list, optional
        A list of catalog definitions from the registry. 
        Defaults to FEUPY_CATALOG_REGISTRY if not provided.

    Returns:
    --------
    List:
        A list of catalog class instances.

    Raises:
    -------
    ValueError: If the catalog class cannot be loaded properly.
    """
    source_catalogs = []
    
    for index, catalog in enumerate(catalogs):
        try:
            catalog_cls = catalogs.get_cls(catalog.tag)()
            source_catalogs.append(catalog_cls)
            # log.info(f"Successfully loaded catalog '{catalog.tag}' at index {index}.")
        except Exception as e:
            log.error(f"Failed to load catalog '{catalog.tag}' at index {index}: {e}")
            raise ValueError(f"Error loading catalog '{catalog.tag}' at index {index}: {e}")
    
    log.info(f"Loaded {len(source_catalogs)} catalogs.")
    return source_catalogs
