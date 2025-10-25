# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Source catalogs."""

from gammapy.catalog import CATALOG_REGISTRY
from gammapy.utils.registry import Registry

from .hawc import (
    get_flux_points_3hwc, create_flux_points_table_3hwc, get_flux_points_2hwc, create_flux_points_table_2hwc,
    SourceCatalogObjectEHWC, SourceCatalogEHWC,
    SourceCatalogObjectExtraHAWC, SourceCatalogExtraHAWC
)
from .hess import (
    SourceCatalogObjectExtraHESS, SourceCatalogExtraHESS
)

from .veritas import SourceCatalogVTSCat, SourceCatalogObjectVTSCat, SourceCatalogVERITAS, SourceCatalogObjectVERITAS
from .psrcat import SourceCatalogPSRCAT, SourceCatalogObjectPSRCAT
from .lhaaso import SourceCatalogObjectLHAASO, SourceCatalogLHAASO
# from .lhaaso import SourceCatalogObject1LHAASO, SourceCatalog1LHAASO
from .lhaaso import SourceCatalogObjectExtraLHAASO, SourceCatalogExtraLHAASO

# Initialize the catalog registry
GAMMAPY_CATALOGS = CATALOG_REGISTRY.copy()

FEUPY_CATALOGS = GAMMAPY_CATALOGS+[
    SourceCatalogEHWC,
    SourceCatalogExtraHAWC,
    SourceCatalogExtraHESS,
    SourceCatalogVTSCat,
    SourceCatalogVERITAS,
    SourceCatalogLHAASO, 
    SourceCatalogExtraLHAASO,
    SourceCatalogPSRCAT, 
]

FEUPY_CATALOG_REGISTRY = Registry(FEUPY_CATALOGS)

__all__ = [
    "FEUPY_CATALOG_REGISTRY",
    "SourceCatalogVTSCat",
    "SourceCatalogObjectVTSCat",
    "SourceCatalogVERITAS",
    "SourceCatalogObjectVERITAS",
    "SourceCatalogPSRCAT",
    "SourceCatalogObjectPSRCAT",
    "SourceCatalogObjectLHAASO",
    "SourceCatalogLHAASO"
    "SourceCatalogExtraLHAASO",
    "SourceCatalogObjectExtraLHAASO",
    "SourceCatalogObjectEHWC",
    "SourceCatalogObjectExtraHAWC",
    "SourceCatalogExtraHAWC",
    "SourceCatalogEHWC",
    "SourceCatalogObjectExtraHESS",
    "SourceCatalogExtraHESS",
    "create_flux_points_table_3hwc",
    "get_flux_points_3hwc",
    "create_flux_points_table_2hwc",
    "get_flux_points_2hwc",
]

