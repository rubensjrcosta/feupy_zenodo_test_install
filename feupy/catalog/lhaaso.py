# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""LHAASO catalog and source classes."""

import logging
import numpy as np
import astropy.units as u
from astropy.table import Table
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling.models import (
    Model, 
    Models, 
    SkyModel, 
    PowerLawSpectralModel, 
    LogParabolaSpectralModel
)
from gammapy.utils.gauss import Gauss2DPDF
from gammapy.utils.scripts import make_path
from gammapy.catalog.core import SourceCatalog, SourceCatalogObject
from gammapy.estimators import FluxPoints
from feupy.utils.table import remove_nan_rows
from feupy.utils.stats import fit_spectral_model_to_flux_points
from feupy.utils.string_handling import string_to_filename_format


# Set up logging
log = logging.getLogger(__name__)

__all__ = [
    "create_flux_points_table_1lhaaso",
    "get_flux_points_1lhaaso",
    "SourceCatalogLHAASO",
    "SourceCatalogObjectLHAASO",
    'SourceCatalogObjectExtraLHAASO',
    'SourceCatalogExtraLHAASO'
]

def create_flux_points_table_1lhaaso(source, which):
    """
    Generate a flux points table for a 1LHAASO catalog source.

    This function extracts differential flux data (`dnde`), reference energy (`e_ref`),
    flux uncertainties (`dnde_err`), and upper limit information from a given source
    in the 1LHAASO catalog. It supports both 'point' and 'extended' source models.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject1LHAASO`
        Source object from the 1LHAASO catalog.
    which : str
        Type of source model to use ('point' or 'extended').
        If 'extended' is specified but the source lacks an extended model, an error is raised.

    Returns
    -------
    table : `~astropy.table.Table`
        Table containing the flux points for the source with the following columns:

        - `e_ref` : Reference energy (TeV).
        - `dnde` : Differential flux at `e_ref`.
        - `dnde_err` : Error on the differential flux.
        - `dnde_ul` : Upper limit for `dnde` (if applicable).
        - `is_ul` : Boolean flag indicating if the data point is an upper limit.

    Raises
    ------
    ValueError
        If 'extended' is specified but the source lacks an extended model.

    Examples
    --------
    >>> from gammapy.catalog import SourceCatalog1LHAASO
    >>> catalog = SourceCatalog1LHAASO()
    >>> source = catalog["1LHAASO J1825-134"]
    >>> table = create_flux_points_table_1lhaaso(source, "point")
    >>> print(table)
    """
    def _parse(source, name, which):
        tag = get_model_tag(source, which)
        is_ul = False
        value = u.Quantity(source.data[f"{name}{tag}"])
        if (
            np.isnan(value) or value == 0 * value.unit
        ) and f"{name}_ul{tag}" in source.data:
            value = source.data[f"{name}_ul{tag}"]
            is_ul = True
        return value, is_ul

    def get_model_tag(source, which):
        if which in source.data["Model_a"]:
            tag = ""
        elif which in source.data["Model_b"]:
            tag = "_b"
        else:
            raise ValueError("Invalid model component name")
        return tag

    def _get(source, name, which):
        value, _ = _parse(source, name, which)
        return value

    data = source.data
    sed_type = 'dnde'
    e_ref = u.Quantity([_get(source, "E0", which)])
    
    spec_model = source.spectral_model(which=which)
    dnde = spec_model(e_ref)
    dnde_err = spec_model.evaluate_error(e_ref)[1]
    
    is_ul = False
    dnde_ul = np.nan
    if not dnde_err.value:
        is_ul = True
        dnde_ul = dnde
        
    table = Table()
    table["e_ref"] = e_ref
    table["e_ref"].description = "Reference energy (TeV)"
    table["dnde"] = dnde
    table["dnde"].description = "Differential flux at reference energy"
    table["dnde_err"] = dnde_err
    table["dnde_err"].description = "Error on the differential flux"
    table["dnde_ul"] = dnde_ul
    table["dnde_ul"].unit = dnde.unit
    table["dnde_ul"].description = "Upper limit for differential flux"
    table["is_ul"] = is_ul
    table["is_ul"].description = "Boolean flag indicating if data is an upper limit"
    
    table.meta['source_name'] = source.name
    table.meta["SED_TYPE"] = sed_type
    table.meta["model"] = which
    table.meta['comments'] = ['Reference: https://iopscience.iop.org/article/10.3847/1538-4365/acfd29']
    
    for column in table.colnames:
        if column.startswith("dnde"):
            table[column].format = ".3e"
        elif column.startswith("e_"):
            table[column].format = ".3f"
    
    return table

def get_flux_points_1lhaaso(source, which):
    """
    Generate flux points as `~gammapy.estimators.FluxPoints` for a 1LHAASO source.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject1LHAASO`
        Source object from the 1LHAASO catalog.
    which : str
        Type of source model ('point' or 'extended').

    Returns
    -------
    flux_points : `~gammapy.estimators.FluxPoints`
        Flux points extracted from the source.
    """
    table = create_flux_points_table_1lhaaso(source, which)
    return FluxPoints.from_table(
        table=table,
        reference_model=source.spectral_model(which),
        sed_type=table.meta['SED_TYPE'],
    )


class SourceCatalogObjectExtraLHAASO(SourceCatalogObject):
    """Represents a single source in the ExtraLHAASO catalog.

    Provides detailed information about a source, including position, spectrum,
    and flux points.

    Attributes
    ----------
    _source_name_key : str
        Name of the key for source names.
    _MODELS : Models
        Pre-loaded models from the extra ExtraLHAASO data file.
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
            "spectrum": self._info_spectrum
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
        models = Models.read("$FEUPY_DATA/data/dedicated_publications/lhaaso/2024icrc.confE.643Y/models.yaml")
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
        file_path = ""
        filename = f"$FEUPY_DATA/data/dedicated_publications/lhaaso/2024icrc.confE.643Y/{string_to_filename_format(self.name)}.fits"
        return FluxPoints.read(filename,  reference_model=self.sky_model(), sed_type='e2dnde')
    
class SourceCatalogExtraLHAASO(SourceCatalog):
    """LHAASO Extra Source Catalog with extended data.

    See: https://doi.org/10.1038/s41586-021-03498-z

    Each source is represented by `SourceCatalogObjectExtraLHAASO`.
    """
    tag = 'LHAASO-2024icrc'
    bibcode = '2024icrc.confE.643Y' 
    description = "LHAASO first 12 PeVatrons Catalogue"

    source_object_class = SourceCatalogObjectExtraLHAASO

    def __init__(self, filename="$FEUPY_DATA/data/dedicated_publications/lhaaso/2024icrc.confE.643Y/catalog.ecsv"):
        table = Table.read(make_path(filename), format='ascii.ecsv')
        super().__init__(table=table, source_name_key="source_name")

class SourceCatalogObjectLHAASO(SourceCatalogObject):
    """One source from the LHAASO first 12 PeVatrons Catalogue.

    See: https://doi.org/10.1038/s41586-021-03498-z

    The data are available through the web page (http://english.ihep.cas.cn/lhaaso/index.html) 
    in the section ‘Public Data’. 

    One source is represented by `~feupy.catalog.SourceCatalogLHAASO`.
    """    
    _source_name_key = "source_name"
    _sed_type = 'e2dnde'
    
    def __str__(self):
        return self.info()
    
    def info(self, info="all"):
        """Summary information string.

        Parameters
        ----------
        info : {'all', 'basic', 'position', 'spectrum'}
            Comma-separated list of options.
        """
        if info == "all":
            info = "basic,position,spectrum"

        ss = ""
        ops = info.split(",")
        if "basic" in ops:
            ss += self._info_basic()
        if "position" in ops:
            ss += self._info_position()
        if "spectrum" in ops:
            ss += self._info_spectrum()

        return ss

    def _info_basic(self):
        """Return basic information about the source."""
        return (
            f"\n*** Basic info ***\n\n"
            f"Catalog row index (zero-based): {self.row_index}\n"
            f"Source name: {self.name}\n"
        )
    
    def _info_position(self):
        """Return position information about the source."""
        return (
            f"\n*** Position info ***\n\n"
            f"RA: {self.data.ra:.3f}\n"
            f"DEC: {self.data.dec:.3f}\n"
            # Uncomment and use if available
            # f"GLON: {self.data.glon:.3f}\n"
            # f"GLAT: {self.data.glat:.3f}\n"
            # f"Position error: {self.data.pos_err:.3f}\n"
        )
    
    def _info_spectrum(self):
        """Return spectral information about the source."""
        ss = "\n*** Spectral info ***\n\n"
        if self.spectral_model is not None:
            model = self.spectral_model()
            parameters = model.parameters
            ss += f"Spectrum type: {model.tag[0]}\n"
            for par in parameters:
                name = par.name
                val = par.value
                err = par.error
                
                try:
                    unit = f"{par.unit:unicode}"
                except AttributeError:
                    unit = ""

                ss += f"{name}: {val:.3f} ± {err} {unit}\n"
        else:
            ss += "No spectrum available"

        return ss
    
    def spectral_model(self):
        """Create and fit a spectral model as a `~gammapy.modeling.models.SpectralModel` object."""
        flux_points_table = self.flux_points_table
        d = self.data
        reference = d['spec_reference']
        spec_type = d['spec_type']

        # Initialize the spectral model based on the type
        if spec_type == 'lp':
            spec_model = LogParabolaSpectralModel(reference=reference)
        elif spec_type == 'pl':
            spec_model = PowerLawSpectralModel(reference=reference)
        else:
            log.warning(f"Unknown spectral model type: {spec_type}")
            return None
        
        return fit_spectral_model_to_flux_points(flux_points_table, spec_model)

    def sky_model(self):
        """Return the source sky model (`~gammapy.modeling.models.SkyModel`)."""
        return SkyModel(
            spectral_model=self.spectral_model(),
            name=self.name,
        )
            
    @property
    def flux_points(self):
        """Return flux points (`~gammapy.estimators.FluxPoints`)."""
        return FluxPoints.from_table(
            table=self.flux_points_table,
            reference_model=self.sky_model(),
            sed_type=self._sed_type,
        )
    
    @property
    def flux_points_table(self):
        """Return differential flux points (`~gammapy.estimators.FluxPoints`)."""
        d = self.data
        spec_type = d['spec_type']
        table = Table()
        table.meta["SED_TYPE"] = self._sed_type
        
        if spec_type == 'pl':
            log.info("Table with a single row generated from the spectral model parameters.")
            
        # Identify valid SED columns that do not contain all NaNs
        valid_sed = [key for key in d.keys() if 'sed' in key and not np.all(np.isnan(d.get(key).value))]
        valid = [key.replace("sed_", "") for key in valid_sed]

        for index, sed_col in enumerate(valid_sed):
            col_name = valid[index]
            if col_name not in table.colnames:
                table[col_name] = d[sed_col]
            else:
                log.warning(f"Column {col_name} already exists in the table.")

        return remove_nan_rows(table)
    

class SourceCatalogLHAASO(SourceCatalog):
    """LHAASO first 12 PeVatrons Catalogue.

    See: https://doi.org/10.1038/s41586-021-03498-z

    The data are available through the web page (http://english.ihep.cas.cn/lhaaso/index.html) 
    in the section ‘Public Data’. 

    One source is represented by `~feupy.catalog.SourceCatalogLHAASO`.
    """    
    tag = "LHAASO"
    bibcode = '2021Natur.594...33C' 
    description = "LHAASO first 12 PeVatrons Catalogue"
    
    source_object_class = SourceCatalogObjectLHAASO
    
    def __init__(self, filename="$FEUPY_DATA/data/catalogs/lhaaso/lhaaso_catalog.ecsv"):
        table = Table.read(make_path(filename), format='ascii.ecsv')
        source_name_key = "source_name"
        super().__init__(table=table, source_name_key=source_name_key)

