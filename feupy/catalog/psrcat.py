# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""ATNF Pulsar Catalogue and source classes."""

import logging
from astropy.table import Table
from gammapy.utils.scripts import make_path
from gammapy.catalog.core import SourceCatalog, SourceCatalogObject

# Set up logging
log = logging.getLogger(__name__)


__all__ = [
    "SourceCatalogPSRCAT",
    "SourceCatalogObjectPSRCAT",
]


class SourceCatalogObjectPSRCAT(SourceCatalogObject):
    """One source from the ATNF Pulsar Catalogue.

    See: Manchester, R. N., Hobbs, G. B., Teoh, A. & Hobbs, M., Astron. J., 129, 1993-2006 (2005) (astro-ph/0412641)

    The data are available through the web page (http://www.atnf.csiro.au/research/pulsar/psrcat)
    in the section ‘Public Data’. 
    """    
    _source_name_key = "NAME"
    
    def __str__(self):
        return self.info()
    
    def info(self, info="all"):
        """Summary information string.

        Parameters
        ----------
        info : {'all', 'basic', 'position', 'timing-profile', 'distance', 'associations-survey', 'derived'}
            Comma-separated list of options for information to include in the summary.
        """
        if info == "all":
            info = "basic,position,timing-profile,distance,associations-survey,derived"

        ss = ""
        ops = info.split(",")
        if "basic" in ops:
            ss += self._info_basic()
        if "position" in ops:
            ss += self._info_position()
        if "timing-profile" in ops:
            ss += self._info_timing_profile()
        if "distance" in ops:
            ss += self._info_distance()
        if "associations-survey" in ops:
            ss += self._info_associations_survey()
        if "derived" in ops:
            ss += self._info_derived()
            
        return ss

    def _info_basic(self):
        """Print basic information."""
        return (
            f"\n*** Basic info ***\n\n"
            f"Catalog row index (zero-based): {self.row_index}\n"
            f"Source name: {self.name}\n"
        )
    
    def _info_position(self):
        """Print position information."""
        return (
            f"\n*** Position info ***\n\n"
            f"RA: {self.data.RAJ2000:.3f} ± {self.data.RAJ2000_ERR:.3f}\n"
            f"DEC: {self.data.DEJ2000:.3f} ± {self.data.DEJ2000_ERR:.3f}\n"
        )
    
    def _info_timing_profile(self):
        """Print timing solution and profile parameters info."""
        return (
            "\n*** Timing and profile info ***\n\n"
            f"P0: {self.data.P0.value:.3e} ± {self.data.P0_ERR:.3e}\n"
        )

    def _info_distance(self):
        """Print distance parameters info."""
        return (
            "\n*** Distance info ***\n\n"
            f"Dist: {self.data.DIST:.2e}\n"
            f"Dist_DM: {self.data.DIST_DM:.2e}\n"
        )
    
    def _info_associations_survey(self):
        """Print associations and survey parameters info."""
        return (
            "\n*** Associations and survey info ***\n\n"
            f"Assoc: {self.data.ASSOC}\n"
            f"Type: {self.data.TYPE}\n"
        )
    
    def _info_derived(self):
        """Print derived parameters info."""
        return (
            "\n*** Derived parameters info ***\n\n"
            f"Age: {self.data.AGE:.2e}\n"
            f"BSurf: {self.data.BSURF:.2e}\n"
            f"E_dot: {self.data.EDOT:.2e}\n"
        )

class SourceCatalogPSRCAT(SourceCatalog):
    """ATNF Pulsar Catalogue.

    See: https://www.atnf.csiro.au/research/pulsar/psrcat/

    One source is represented by `SourceCatalogObjectPSRCAT`.
    """  
    
    tag = "psrcat"
    bibcode = '2005AJ....129.1993M'    
    description = "ATNF Pulsar Catalogue, a comprehensive database of all published pulsars"
        
    source_object_class = SourceCatalogObjectPSRCAT
    
    def __init__(self, filename="$FEUPY_DATA/data/catalogs/psrcat/psrcat_catalog.fits"):
        try:
            table = Table.read(make_path(filename))
        except Exception as e:
            log.error(f"Error loading the catalog: {e}")
            raise

        super().__init__(table=table, source_name_key="NAME")
    
    @property
    def PSR_PARAMS(self):
        return self.table.colnames
    
    @property
    def PSR_PARAMS_DESCRIPTION(self):
        """Returns the description of pulsar parameters."""
        ss = "\n*** The Pulsar Parameters ***\n\n" 
        for par in self.PSR_PARAMS:
            unit = f" ({self.table[par].unit})" if self.table[par].unit is not None else ""
            ss += f"{self.table[par].name}: {self.table[par].description}{unit}\n"
        return ss
