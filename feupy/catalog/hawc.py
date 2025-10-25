# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""HAWC Source Catalog."""

import numpy as np
from astropy import units as u
from astropy.table import Table, Column
from gammapy.utils.scripts import make_path
from gammapy.estimators import FluxPoints
from gammapy.modeling.models import SkyModel, Models
from gammapy.catalog.core import SourceCatalog, SourceCatalogObject
from gammapy.catalog.hawc import SourceCatalog3HWC, SourceCatalog2HWC
import logging
from feupy.utils.string_handling import string_to_filename_format

# Set up logging
log = logging.getLogger(__name__)

__all__ = [
    "create_flux_points_table_3hwc",
    "get_flux_points_3hwc",
    "create_flux_points_table_2hwc",
    "get_flux_points_2hwc",
    "SourceCatalogObjectEHWC",
    "SourceCatalogEHWC",
    "SourceCatalogObjectExtraHAWC",
    "SourceCatalogExtraHAWC",
    
]
        
def create_flux_points_table_3hwc(source):
    """
    Generate a flux points table for a 3HWC catalog source.

    This function creates a flux points table containing differential flux
    data (`dnde`), reference energy (`e_ref`), flux uncertainties (`dnde_errn`, `dnde_errp`),
    and upper limit information for a specified source in the 3HWC catalog.
    
    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject3HWC`
        Source object from the 3HWC catalog for which the flux points table
        is to be created.
        
    Returns
    -------
    table : `~astropy.table.Table`
        Table containing flux points for the source, with metadata and columns:
        
        - `e_ref` : Reference energy for the differential flux data (TeV).
        - `dnde` : Differential flux (`dnde`) at `e_ref`.
        - `dnde_errn` : Negative error on `dnde`.
        - `dnde_errp` : Positive error on `dnde`.
        - `dnde_ul` : Upper limit for `dnde` (if applicable).
        - `is_ul` : Boolean flag indicating if the data point is an upper limit.

    Notes
    -----
    The `dnde_ul` column contains NaN if no upper limit is provided for the source,
    and `is_ul` is set to `False` by default.

    Examples
    --------
    >>> from gammapy.catalog import SourceCatalog3HWC
    >>> catalog = SourceCatalog3HWC()
    >>> source = catalog["3HWC J1825-134"]
    >>> table = create_flux_points_table_3hwc(source)
    >>> print(table)
    """
    # Initialize catalog and source data
    catalog = SourceCatalog3HWC()
    data = source.data
    sed_type = "dnde"
    
    # Create a table and add metadata
    table = Table()
    table.meta['source_name'] = source.name
    table.meta['catalog_name'] = catalog.table.meta['catalog_name']
    table.meta["SED_TYPE"] = sed_type
    table.meta['search_radius'] = data['search_radius']
    table.meta['spec0_radius'] = data['spec0_radius']
    table.meta['reference'] = catalog.table.meta['reference']
    
    # Define and add columns
    col_1 = Column(
        name='e_ref',
        data=u.Quantity([7*u.TeV]),
        description='Reference Energy',
        format='.3g'
    )
    
    col_2 = Column(
        name='dnde',
        data=u.Quantity([data['spec0_dnde']]),
        description=catalog.table['spec0_dnde'].description,
        format=catalog.table['spec0_dnde'].format
    )
    
    col_3 = Column(
        name='dnde_errn',
        data=u.Quantity([data['spec0_dnde_errn'] * -1]),
        description=catalog.table['spec0_dnde_errn'].description,
        format=catalog.table['spec0_dnde_errn'].format
    )
    
    col_4 = Column(
        name='dnde_errp',
        data=u.Quantity([data['spec0_dnde_errp']]),
        description=catalog.table['spec0_dnde_errp'].description,
        format=catalog.table['spec0_dnde_errp'].format
    )
    
    col_5 = Column(
        name='dnde_ul',
        data=[np.nan],
        description='Differential flux (dnde) SED upper limit',
        unit=data['spec0_dnde'].unit
    )
    
    col_6 = Column(
        name='is_ul',
        data=[False],
        description='Whether data is an upper limit.',
        dtype=bool
    )
    
    # Add columns to the table
    table.add_columns([col_1, col_2, col_3, col_4, col_5, col_6])
    
    return table


def get_flux_points_3hwc(source):
    """
    Generate `FluxPoints` for a 3HWC catalog source.

    Uses the flux points table for a 3HWC source to create a `FluxPoints`
    object with associated spectral model.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject3HWC`
        Source object from the 3HWC catalog.

    Returns
    -------
    flux_points : `~gammapy.estimators.FluxPoints`
        Flux points object containing the flux data for the source.
    
    Examples
    --------
    >>> from gammapy.catalog import SourceCatalog3HWC
    >>> catalog = SourceCatalog3HWC()
    >>> source = catalog["3HWC J1825-134"]
    >>> flux_points = get_flux_points_3hwc(source)
    >>> print(flux_points)
    """
    table = create_flux_points_table_3hwc(source)
    spec_model = source.spectral_model()
    return FluxPoints.from_table(
        table, 
        sed_type=table.meta["SED_TYPE"],
        reference_model=spec_model
    )


def create_flux_points_table_2hwc(source, which='point'):
    """
    Generate a flux points table for a 2HWC catalog source.

    Creates a flux points table containing differential flux data (`dnde`),
    reference energy (`e_ref`), flux uncertainties (`dnde_err`), and upper
    limit information for a specified source in the 2HWC catalog. The function
    handles both 'point' and 'extended' models.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject2HWC`
        Source object from the 2HWC catalog for which the flux points table
        is to be created.
    which : str, optional
        Type of source model to use ('point' or 'extended'). Default is 'point'.
        If 'extended' is specified and the source lacks an extended model,
        an error is raised.

    Returns
    -------
    table : `~astropy.table.Table`
        Table containing the flux points for the given source with columns:
        
        - `e_ref` : Reference energy (TeV).
        - `dnde` : Differential flux at `e_ref`.
        - `dnde_err` : Error on the differential flux.
        - `dnde_ul` : Upper limit for `dnde` (if applicable).
        - `is_ul` : Boolean flag indicating if the data point is an upper limit.

    Raises
    ------
    ValueError
        If 'extended' is specified for `which` and the source lacks an extended model.

    Examples
    --------
    >>> from gammapy.catalog import SourceCatalog2HWC
    >>> catalog = SourceCatalog2HWC()
    >>> source = catalog["2HWC J1825-134"]
    >>> table = create_flux_points_table_2hwc(source)
    >>> print(table)
    """
    # Initialize catalog and source data
    catalog = SourceCatalog2HWC()
    data = source.data
    sed_type = "dnde"
    
    # Check for extended model if specified
    if which == 'extended' and source.n_models != 2:
        raise ValueError("No extended model available for this source.")
    
    # Define reference energy and model evaluation
    e_ref = u.Quantity([7 * u.TeV])
    spec_model = source.spectral_model(which=which)
    dnde = spec_model(e_ref)
    dnde_err = spec_model.evaluate_error(e_ref)[1]
    
    # Determine upper limit
    is_ul = False
    dnde_ul = np.nan
    if not dnde_err.value:
        is_ul = True
        dnde_ul = dnde

    # Create the flux points table and add data
    table = Table()
    table["e_ref"] = e_ref
    table["e_ref"].description = 'Reference energy'
    table["dnde"] = dnde
    table["dnde"].description = 'Differential flux (dnde) SED values'
    table["dnde_err"] = dnde_err
    table["dnde_err"].description = 'Differential flux (dnde) SED errors'
    table["dnde_ul"] = dnde_ul
    table["dnde_ul"].description = 'Differential flux (dnde) SED upper limit'
    table["is_ul"] = is_ul
    table["is_ul"].description = 'Indicates if data is an upper limit.'

    # Add metadata to the table
    table.meta['source_name'] = source.name
    table.meta['catalog_name'] = catalog.table.meta['catalog_name']
    table.meta["SED_TYPE"] = sed_type
    table.meta['reference'] = catalog.table.meta['reference']

    # Format columns for readability
    for column in table.colnames:
        if column.startswith("dnde"):
            table[column].format = ".3e"
        elif column.startswith("e_"):
            table[column].format = ".3f"

    return table

def get_flux_points_2hwc(source, which='point'):
    """
    Generate `FluxPoints` for a 2HWC catalog source.

    This function creates a `FluxPoints` object for a source in the 2HWC catalog, 
    using the flux points table and associated spectral model. The function can handle 
    both 'point' and 'extended' models.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject2HWC`
        Source object from the 2HWC catalog.
    which : str, optional
        Type of source model to use ('point' or 'extended'). Default is 'point'.
        If 'extended' is specified and the source lacks an extended model,
        an error is raised.

    Returns
    -------
    flux_points : `~gammapy.estimators.FluxPoints`
        Flux points object containing the flux data for the source, including 
        differential flux values and associated uncertainties.

    Raises
    ------
    ValueError
        If 'extended' is specified for `which` and the source lacks an extended model.

    Examples
    --------
    >>> from gammapy.catalog import SourceCatalog2HWC
    >>> catalog = SourceCatalog2HWC()
    >>> source = catalog["2HWC J1825-134"]
    >>> flux_points = get_flux_points_2hwc(source)
    >>> print(flux_points)
    """
    # Generate flux points table for the specified source and model type
    table = create_flux_points_table_2hwc(source, which=which)
    
    # Retrieve the spectral model for the source
    spec_model = source.spectral_model(which=which)
    
    # Create and return the FluxPoints object
    return FluxPoints.from_table(
        table, 
        sed_type=table.meta["SED_TYPE"],
        reference_model=spec_model
    )

class SourceCatalogObjectEHWC(SourceCatalogObject):
    """Represents a single source in the HAWC catalog.

    Provides detailed information about a source, including position, spectrum,
    and flux points.

    Attributes
    ----------
    _source_name_key : str
        Name of the key for source names.
    _MODELS : Models
        Pre-loaded models from the extra HAWC data file.
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
        models = Models.read("$FEUPY_DATA/data/catalogs/ehwc/models.yaml")
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
        if self.flux_points_table:
            return FluxPoints.from_table(
                table=self.flux_points_table,
    #             reference_model=self.sky_model(),
                sed_type='e2dnde',
            )
        return
    
    def _add_source_meta(self, table):
        """Copy over some information to `table.meta`."""
        d = self.data
        m = table.meta
        catalog = SourceCatalogEHWC()
        m["source_name"] = self.name
        m["catalog_name"] = catalog.table.meta['catalog_name']
        m["SED_TYPE"] = "e2dnde"
        m["comments"] = catalog.table.meta['comments']        
        
    @property
    def flux_points_table(self):
        """Differential flux points (`~gammapy.estimators.FluxPoints`)."""
        
        d = self.data
        
        table = Table()
        
        # Set metadata for the table
        self._add_source_meta(table)
        
        valid = np.isfinite(d["sed_e_ref"].value)

        if valid.sum() == 0:
            return None

        table["e_ref"] = d["sed_e_ref"]
        table["e_ref"].description = 'Reference energy'
        table["e2dnde"] = d["sed_e2dnde"]
        table["e2dnde"].description = 'Differential flux (e2dnde) SED values'
        table["e2dnde_errn"] = d["sed_e2dnde_errn"]
        table["e2dnde_errn"].description = 'Differential flux (e2dnde) SED negative errors'
        table["e2dnde_errp"] = d["sed_e2dnde_errp"]
        table["e2dnde_errp"].description = 'Differential flux (e2dnde) SED positive errors'
        table["e2dnde_ul"] = d["sed_e2dnde_ul"]
        table["e2dnde_ul"].description = 'Differential flux (e2dnde) SED upper limit'
        table["is_ul"] = d["sed_is_ul"]        
        table["is_ul"].description = 'Upper limit indicator'

        # Format numeric columns
        for col in table.colnames:
            if col.startswith("e2dnde"):
                table[col].format = ".3e"
            elif col.startswith("e_"):
                table[col].format = ".3f"
                
        # Only keep rows that actually contain information
        table = table[valid]

        for colname in table.colnames:
            if not np.isfinite(table[colname]).any():
                table.remove_column(colname)
        return table


class SourceCatalogEHWC(SourceCatalog):
    """HAWC Extra Source Catalog with extended data.

    See: https://doi.org/10.1103/PhysRevLett.124.021102

    Each source is represented by `SourceCatalogObjectHAWC`.
    """
    tag = "ehwc"
    bibcode = '2020PhRvL.124b1102A' 
    description = "Extra HAWC catalog data"

    source_object_class = SourceCatalogObjectEHWC

    def __init__(self, filename="$FEUPY_DATA/data/catalogs/ehwc/ehwc_catalog.ecsv"):
        table = Table.read(make_path(filename), format='ascii.ecsv')
        super().__init__(table=table, source_name_key="source_name")
        
class SourceCatalogObjectExtraHAWC(SourceCatalogObject):
    """Represents a single source in the ExtraHAWC catalog.

    Provides detailed information about a source, including position, spectrum,
    and flux points.

    Attributes
    ----------
    _source_name_key : str
        Name of the key for source names.
    _MODELS : Models
        Pre-loaded models from the extra ExtraHAWC data file.
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
        models = Models.read("$FEUPY_DATA/data/dedicated_publications/hawc/2021ApJ...907L..30A/models.yaml")
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
        filename = f"$FEUPY_DATA/data/dedicated_publications/hawc/2021ApJ...907L..30A/{string_to_filename_format(self.name)}.fits"
        return FluxPoints.read(filename,  reference_model=self.sky_model(), sed_type='e2dnde')
    
class SourceCatalogExtraHAWC(SourceCatalog):
    """LHAASO Extra Source Catalog with extended data.

    See: https://iopscience.iop.org/article/10.3847/2041-8213/abd77b

    Each source is represented by `SourceCatalogObjectExtraHAWC`.
    """
    tag = 'hwc-2021ApJ'
    bibcode = '2021ApJ...907L..30A' 
    description = "Evidence of 200 TeV Photons from HAWC J1825-134 "

    source_object_class = SourceCatalogObjectExtraHAWC

    def __init__(self, filename="$FEUPY_DATA/data/dedicated_publications/hawc/2021ApJ...907L..30A/catalog.ecsv"):
        table = Table.read(make_path(filename), format='ascii.ecsv')
        super().__init__(table=table, source_name_key="source_name")
