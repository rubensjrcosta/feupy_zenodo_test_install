# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""VTSCat and source classes."""

import os
import logging
import string

import numpy as np
import pandas as pd
from pandas import json_normalize
from astropy.table import Table

from gammapy.modeling.models import (Model, Models, SkyModel, 
                                     PowerLawSpectralModel, LogParabolaSpectralModel)
from gammapy.modeling import Fit
from gammapy.datasets import Datasets, FluxPointsDataset
from gammapy.estimators import FluxPoints
from gammapy.catalog.core import SourceCatalog, SourceCatalogObject
from gammapy.utils.scripts import make_path

from feupy.utils.io import read_yaml
from feupy.utils.string_handling import string_to_filename_format

# Set up logging
log = logging.getLogger(__name__)

__all__ = [
    "SourceCatalogVTSCat",
    "SourceCatalogObjectVTSCat",
    "SourceCatalogVERITAS",
    "SourceCatalogObjectVERITAS",
]

def generate_unique_name(name, reference_id, unique_names):
    """Generate a unique name by appending reference_id and, if necessary, a letter."""
    new_name = f'{name} ({reference_id[0:7]})'
    if new_name in unique_names:
        for letter in string.ascii_letters:
            new_name = f'{name} ({reference_id[0:7]}-{letter})'
            if new_name not in unique_names:
                return new_name
    else:
        return new_name

class SourceCatalogObjectVTSCat(SourceCatalogObject):
    """One source from the VTSCat Catalogue.

    The data are available through the web page (https://iopscience.iop.org/article/10.3847/1538-4357/aac4a2) 
    in tables 8-13.

    One source is represented by `~feupy.catalog.SourceCatalogVTSCat`.
    """    

    _DATASETS_PATH = '$FEUPY_DATA/data/catalogs/vtscat/datasets'
    _source_name_key = "source_name"
        
    def __str__(self):
        return self.info()
    
    def info(self, info="all"):
        """Summary information string.

        Parameters
        ----------
        info : {'all', 'basic', 'position', 'spectrum'}
            Comma-separated list of options
        """
        if info == "all":
            info = "basic,position,spectrum"

        ss = ""
        ops = info.split(",")
        if "basic" in ops:
            ss += self._info_basic()
        if "position" in ops:
            ss += self._info_position()

        return ss

    def _info_basic(self):
        """Return basic information about the source."""
        d = self.data
        return (
            f"\n*** Basic info ***\n\n"
            f"Catalog row index (zero-based): {self.row_index}\n"
            f"Source name: {self.name}\n"
            f"VTSCat name: {d.veritas_name}\n"
            f"Common name: {d.common_name}\n"
            f"Other names: {d.other_names}\n"
            f"VTSCat id: {d.veritas_id}\n"
            f"Location: {d.where}\n"
            f"Type: {d.type}\n"
            f"VTSCat components: {d.veritas_components}\n"
            f"VTSCat ID name: {self.veritas_id}\n"
            f"Simbad ID: {d.simbad_id}\n"
            f"References: {d.reference_id}\n\n"
        )
    
    def _info_position(self):
        """Return position information about the source."""
        return (
            f"\n*** Position info ***\n\n"
            f"RA: {self.data.ra:.3f}\n"
            f"DEC: {self.data.dec:.3f}\n"
        )

    def _spectral_model(self, table):
        """Create and fit a spectral model as a `~gammapy.modeling.models.SpectralModel` object."""
        reference = "1 TeV"
        spec_type = 'pl'

        if spec_type == 'lp':
            spec_model = LogParabolaSpectralModel(reference=reference)
        elif spec_type == 'pl':
            spec_model = PowerLawSpectralModel(reference=reference)
        else:
            log.warning(f"Unknown spectral model type: {spec_type}")
            return None
            
        flux_points = FluxPoints.from_table(table)
        ds = FluxPointsDataset(data=flux_points)
        datasets = Datasets(ds)

        model = SkyModel(spectral_model=spec_model)
        datasets.models = model
        fitter = Fit()
        result = fitter.run(datasets=datasets)
        
        return model.spectral_model

    def _sky_model(self, table):
        """Return the source sky model (`~gammapy.modeling.models.SkyModel`)."""
        _models_names = []
        meta = table.meta
        reference_id = meta.get('reference_id')
        
        model_name = generate_unique_name(self.name, reference_id, _models_names)
        _models_names.append(model_name)
        
        return SkyModel(
            spectral_model=self._spectral_model(table),
            name=model_name,
        )
            
    def flux_points(self):
        """Return flux points tables as a list of Astropy Tables."""
        tables = []
        for table in self.flux_points_tables:
            tables.append(self._flux_points(table))  
        return tables
    
    def _flux_points(self, table):
        """Return flux points (`~gammapy.estimators.FluxPoints`)."""
#         print(table.meta)
        return FluxPoints.from_table(
            table=table,
            reference_model=self._sky_model(table),
            sed_type=table.meta.get('SED_TYPE'),
        )

    def _reference_id(self):
        """Return reference IDs as a cleaned list."""
        _ids = self.data['reference_id'].split(", ")
        return [_.replace(" ", "") for _ in _ids]
    
    @property
    def veritas_id(self):
        """Return the formatted VERITAS ID."""
        return f'VER-{self.data.veritas_id:06}'
        
    def _get_file_paths(self, reference_id, which='info'):
        """Get file paths based on the reference IDs."""
        paths = []
        yaml_dir = make_path(f'{self._DATASETS_PATH}/{reference_id}')
        for filename in os.listdir(yaml_dir):
            path = make_path(f'{yaml_dir}/{filename}')
            if which == 'info' and filename == 'info.yaml':
                paths.append(path)
            elif which in ['sed', 'lc']:
                if filename.endswith('.ecsv') and self.veritas_id in filename and which in filename:
                    paths.append(path)
            elif which == 'obs' and filename.endswith('.yaml') and self.veritas_id in filename:
                paths.append(path)
        return paths
    
    @property
    def reference_id_info(self):
        """Return reference IDs as a normalized JSON object."""
        data = []
        for ref in self._reference_id():
            for path in self._get_file_paths(ref, which='info'):
                data.append(read_yaml(path))
        return json_normalize(data)

    def get_observation_tables(self, reference_id):
        """Get file paths based on the reference IDs."""
        data = []        
        for path in self._get_file_paths(reference_id, which='obs'):                
            data.append(read_yaml(path))                  
        return data
        
    @property
    def observation_info(self):
        """Return observation tables as a normalized JSON object."""
        data = []
        for ref in self._reference_id():
            data.extend(self.get_observation_tables(ref))                
        return json_normalize(data)

    def get_flux_points_tables(self, reference_id):
        """Get file paths based on the reference IDs."""
        tables = []        
        for path in self._get_file_paths(reference_id, which='sed'):                
            table = Table.read(make_path(path), format='ascii.ecsv')
            if 'dnde' in table.colnames:
                table.meta["SED_TYPE"] = 'dnde'
                if 'dnde_ul' in table.colnames:
                    table['is_ul'] = [not np.isnan(_) for _ in table['dnde_ul']]
            if 'e2dnde' in table.colnames:
                table.meta["SED_TYPE"] = 'e2dnde'
                if 'e2dnde_ul' in table.colnames:
                    table['is_ul'] = [not np.isnan(_) for _ in table['e2dnde_ul']]  
            tables.append(table)
        return tables
    
    @property
    def flux_points_tables(self):
        """Return flux points tables as a list of Astropy Tables."""
        tables = []
        for ref in self._reference_id():
            tables.extend(self.get_flux_points_tables(ref))  
        return tables

class SourceCatalogVTSCat(SourceCatalog):
    """VTSCat Catalogue.
        One source is represented by `~feupy.catalog.SourceCatalogVTSCat`.
        See https://iopscience.iop.org/article/10.3847/2515-5172/acb147
    """    
    tag = "vtscat"
    bibcode = '2023RNAAS...7....6A' 
    description = "VTSCat catalog from the VTSCat observatory"
    
    source_object_class = SourceCatalogObjectVTSCat
    
    def __init__(self, filename="$FEUPY_DATA/data/catalogs/vtscat/sources/vtscat.ecsv"):
        table = Table.read(make_path(filename), format='ascii.ecsv')
        source_name_key = "source_name"
        source_name_alias = ("veritas_name", "common_name", "other_names", "simbad_id")
        super().__init__(table=table, source_name_key=source_name_key)

class SourceCatalogObjectVERITAS(SourceCatalogObject):
    """One source from the VERITAS Catalogue.    
    The data are available through the web page (https://iopscience.iop.org/article/10.3847/1538-4357/aac4a2) 
    in the tables 8-13. 

    One source is represented by `~feupy.catalog.SourceCatalogVERITAS`.
    """    
    _source_name_key = "source_name"
    
    _DATA_PATH = "$FEUPY_DATA/data/catalogs/veritas/"  
    _MODELS = Models.read(f"{_DATA_PATH}/models.yaml")   

    def __str__(self):
        return self.info()
    
    def info(self, info="all"):
        """Summary information string.

        Parameters
        ----------
        info : {'all', 'basic', 'position', 'spectrum'}
            Comma separated list of options
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
        """Print basic information."""
        return (
            f"\n*** Basic info ***\n\n"
            f"Catalog row index (zero-based) : {self.row_index}\n"
            f"Source name : {self.name}\n"
        )
    
    def _info_position(self):
        """Print position information."""
        return (
            f"\n*** Position info ***\n\n"
            f"RA: {self.data.ra:.3f}\n"
            f"DEC: {self.data.dec:.3f}\n"
        )
    
    def _info_spectrum(self):
        """Print spectral info."""
        ss = "\n*** Spectral info ***\n\n"
        if self.spectral_model is not None:
            model = self.spectral_model()
            parameters = model.parameters
            ss += f"Spectrum type:  {model.tag[0]}\n"
            for par in parameters:
                name = par.name
                val = par.value
                err = par.error
                
                try:
                    unit = f"{par.unit:unicode}"
                except: unit = ""

                ss += f"{name}: {val:.3} +- {err:.3} {unit}\n"
    
        else:
            ss += "No spectrum available"

        return ss
  
    def spectral_model(self):
        """Spectral model as a `~gammapy.modeling.models.SpectralModel` object."""
        return self._MODELS[self.name].spectral_model
    
    def spatial_model(self):
        """Spatial model as a `~gammapy.modeling.models.SpatialModel` object."""
        return self._MODELS[self.name].spatial_model
    
    
    def sky_model(self):
        """Source sky model (`~gammapy.modeling.models.SkyModel`)."""
        return self._MODELS[self.name] 
    
    @property
    def flux_points(self):
        """Flux points (`~gammapy.estimators.FluxPoints`)."""
        filename = f'{self._DATA_PATH}/{string_to_filename_format(self.name)}.fits'
        return FluxPoints.read(filename)
    
class SourceCatalogVERITAS(SourceCatalog):
    """VERITAS  Catalogue.

    See: https://iopscience.iop.org/article/10.3847/1538-4357/aac4a2

    One source is represented by `~feupy.catalog.SourceCatalogVERITAS`.
    """    
    tag = "veritas-2018ApJ"
    bibcode = '2018ApJ...861..134A'         
    description = "A Very High Energy γ-Ray Survey toward the Cygnus Region of the Galaxy"
    
    source_object_class = SourceCatalogObjectVERITAS
    
    def __init__(self, filename="$FEUPY_DATA/data/catalogs/veritas/veritas.fits"):
        table = Table.read(make_path(filename))
        source_name_key = "source_name"
        super().__init__(table=table, source_name_key=source_name_key)
        
