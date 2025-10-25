# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""HESS Source Catalog."""

from astropy.table import Table
from gammapy.utils.scripts import make_path
from gammapy.estimators import FluxPoints
from gammapy.modeling.models import SkyModel, Models
from gammapy.catalog.core import SourceCatalog, SourceCatalogObject
from feupy.utils.string_handling import string_to_filename_format
import logging

# Set up logging
log = logging.getLogger(__name__)


__all__ = [
    "SourceCatalogObjectExtraHESS",
    "SourceCatalogExtraHESS",
]


class SourceCatalogObjectExtraHESS(SourceCatalogObject):
    """Represents a single source in the ExtraHESS catalog.

    Provides detailed information about a source, including position, spectrum,
    and flux points.

    Attributes
    ----------
    _source_name_key : str
        Name of the key for source names.
    _MODELS : Models
        Pre-loaded models from the ExtraHESS data file.
    """
    _source_name_key = "source_name"
    
    def __str__(self):
        return self.info()

    def info(self, info="all"):
        """Summary information string for the source.

        Parameters
        ----------
        info : {'all', 'basic', 'position', 'spectrum'}
            Type of information to display. Options are:
            - 'all': Shows basic, position, and spectrum information.
            - 'basic': Shows basic source information.
            - 'position': Shows position information.
            - 'spectrum': Shows spectral information.
        """
        details = {
            "basic": self._info_basic,
            "position": self._info_position,
            "spectrum": self._info_spectrum,
        }
        selected_info = info.split(",") if info != "all" else details.keys()
        return "\n".join(details[opt]() for opt in selected_info if opt in details)

    def _info_basic(self):
        return f"\n*** Basic info ***\n\nCatalog row index: {self.row_index}\nSource name: {self.name}\n"

    def _info_position(self):
        return f"\n*** Position info ***\n\nRA: {self.data.ra:.3f}\nDEC: {self.data.dec:.3f}\n"

    def _info_spectrum(self):
        if self.spectral_model() is None:
            return "No spectral information available."

        model = self.spectral_model()
        return "\n".join([
            "\n*** Spectral info ***\n",
            f"Spectrum type: {model.tag[0]}",
            *(f"{par.name}: {par.value:.3f} ± {par.error} {par.unit if par.unit else ''}" for par in model.parameters)
        ])

    def spectral_model(self):
        """Get the spectral model associated with this source."""
        models = Models.read("$FEUPY_DATA/data/dedicated_publications/hess/2019Apercent26A...621A.116H/models.yaml")
        if self.name in models.names:
            return models[self.name].spectral_model
        return None

    def sky_model(self):
        """Create a SkyModel representation of the source."""
        if self.spectral_model():
            return SkyModel(spectral_model=self.spectral_model(), name=self.name)
        return None

    @property
    def flux_points(self):
        """Flux points as a `~gammapy.estimators.FluxPoints` object."""
        filename = f"$FEUPY_DATA/data/dedicated_publications/hess/2019Apercent26A...621A.116H/{string_to_filename_format(self.name)}.fits"
        return FluxPoints.read(filename, reference_model=self.sky_model(), sed_type='e2dnde')
    
class SourceCatalogExtraHESS(SourceCatalog):
    """HESS Extra Source Catalog with extended data.

    See: https://www.aanda.org/articles/aa/full_html/2019/01/aa34335-18/aa34335-18.html

    Each source is represented by `SourceCatalogObjectExtraHESS`.
    """
    tag = 'hess-2019A&A'
    bibcode = '2019A&A...621A.116H'
    description = "Particle transport within the pulsar wind nebula HESS J1825–137★"

    source_object_class = SourceCatalogObjectExtraHESS

    def __init__(self, filename="$FEUPY_DATA/data/dedicated_publications/hess/2019Apercent26A...621A.116H/catalog.ecsv"):
        table = Table.read(make_path(filename), format='ascii.ecsv')
        super().__init__(table=table, source_name_key="source_name")
