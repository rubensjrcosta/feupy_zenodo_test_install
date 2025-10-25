# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities for analysis configuration."""

from typing import List, Optional
from gammapy.analysis.config import (
GammapyBaseConfig, SpatialCircleConfig, GeomConfig, BackgroundConfig, SafeMaskConfig, FitConfig,
)
from gammapy.makers import MapDatasetMaker
from gammapy.utils.types import QuantityType, AngleType, PathType
from feupy.utils.types import IrfType
from feupy.utils.enum import TableEnum
from gammapy.analysis.config import (
    MapSelectionEnum,
    ReductionTypeEnum,
)


__all__ = [
    "ObservationConfig",
    "DatasetsConfig",
    "OnOffConfig",
    "SensitivityConfig",
    "StatisticsConfig",
]

class OnOffConfig(GammapyBaseConfig):
    acceptance: int = 1
    acceptance_off: int = 5
    
class DatasetsConfig(GammapyBaseConfig):
    type: ReductionTypeEnum = ReductionTypeEnum.spectrum
    stack: bool = True
    geom: GeomConfig = GeomConfig()
    map_selection: List[MapSelectionEnum] = MapDatasetMaker.available_selection
    background: BackgroundConfig = BackgroundConfig()
    safe_mask: SafeMaskConfig = SafeMaskConfig()
    on_region: SpatialCircleConfig = SpatialCircleConfig()
    containment_correction: bool = True
    use_region_center: bool = False
    on_off: OnOffConfig = OnOffConfig()
    containment: float = 0.68

class StatisticsConfig(GammapyBaseConfig):
    # alpha: float = None
    n_obs: int = 1
    # wstat: dict = {}
    # fitted_parameters: dict = {}  

class SensitivityConfig(GammapyBaseConfig):
    gamma_min: int = 10
    n_sigma: int = 5
    bkg_syst_fraction: float = 0.05
    data_path: Optional[PathType] = None
    table_format: TableEnum = 'fits'

class ObservationConfig(GammapyBaseConfig):
    obs_cone: SpatialCircleConfig = SpatialCircleConfig()
    livetime: Optional[QuantityType] = None
    offset: Optional[QuantityType] = None
    position_angle: Optional[AngleType] = None
    required_irfs: IrfType = ['South', 'AverageAz', '20deg', '50h']
    
# from pydantic.v1 import BaseModel
# from astropy.coordinates import Angle
# from astropy.units import Quantity
# from astropy.time import Time
# from typing import List
# from pathlib import Path

# from feupy.utils.types import (
#     AngleType,
#     EnergyType,
#     QuantityType,
#     TimeType,
#     IrfType,
# )
# from feupy.utils.enum import (
#     FrameEnum,
#     BackgroundMethodEnum,
#     TableEnum,
#     SafeMaskMethodsEnum,
#     MapSelectionEnum,
#     RequiredHDUEnum,
#     ReductionTypeEnum,
# )
# from gammapy.makers import MapDatasetMaker


# # Base Configuration
# class FeupyBaseConfig(BaseModel):
#     class Config:
#         validate_all = True
#         validate_assignment = True
#         arbitrary_types_allowed = True
#         from_attributes = True
#         extra = "forbid"
#         json_encoders = {
#             Angle: lambda v: f"{v.value} {v.unit}",
#             Quantity: lambda v: f"{v.value} {v.unit}",
#             Time: lambda v: f"{v.value}",
#         }


# # Logging Configuration
# class LogConfig(FeupyBaseConfig):
#     level: str = "info"
#     filename: Path = None
#     filemode: str = None
#     format: str = None
#     datefmt: str = None


# # General Configuration
# class GeneralConfig(FeupyBaseConfig):
#     log: LogConfig = LogConfig()
#     outdir: str = "."
#     n_jobs: int = 1
#     data_path: Path = None
#     config_file: Path = None
#     datasets_file: Path = None
#     models_file: Path = None
#     sources_file: Path = None
    

# # Sky and Positioning Configurations
# class SkyCoordConfig(FeupyBaseConfig):
#     frame: FrameEnum = None
#     lon: AngleType = None
#     lat: AngleType = None


# class ROIConfig(FeupyBaseConfig):
#     name: str = 'ROI'
#     position: SkyCoordConfig = SkyCoordConfig()
#     radius: AngleType = None


# class TargetConfig(FeupyBaseConfig):
#     name: str = None
#     position: SkyCoordConfig = SkyCoordConfig()
#     model: dict = {}
#     redshift: float = None


# # Energy Configuration
# class EnergyRangeConfig(FeupyBaseConfig):
#     min: EnergyType = None
#     max: EnergyType = None


# class EnergyAxisConfig(FeupyBaseConfig):
#     min: EnergyType = None
#     max: EnergyType = None
#     nbins: int = None
#     name: str = "energy"


# class EnergyAxisTrueConfig(FeupyBaseConfig):
#     min: EnergyType = None
#     max: EnergyType = None
#     nbins: int = None
#     name: str = "energy_true"


# class EnergyAxesConfig(FeupyBaseConfig):
#     energy: EnergyAxisConfig = EnergyAxisConfig()
#     energy_true: EnergyAxisTrueConfig = EnergyAxisTrueConfig()


# # Observation Configuration
# class ObservationConfig(FeupyBaseConfig):
#     target: TargetConfig = TargetConfig()
#     livetime: QuantityType = None
#     offset: QuantityType = None
#     position_angle: AngleType = None
#     required_irfs: IrfType = ['South', 'AverageAz', '20deg', '50h']


# # Geometry and Map Configuration
# class WidthConfig(FeupyBaseConfig):
#     width: AngleType = "5 deg"
#     height: AngleType = "5 deg"


# class WcsConfig(FeupyBaseConfig):
#     skydir: SkyCoordConfig = SkyCoordConfig()
#     binsize: AngleType = "0.02 deg"
#     width: WidthConfig = WidthConfig()
#     binsize_irf: AngleType = "0.2 deg"


# class SelectionConfig(FeupyBaseConfig):
#     offset_max: AngleType = "2.5 deg"


# class GeomConfig(FeupyBaseConfig):
#     wcs: WcsConfig = WcsConfig()
#     selection: SelectionConfig = SelectionConfig()
#     axes: EnergyAxesConfig = EnergyAxesConfig()

# class FitConfig(FeupyBaseConfig):
#     fit_range: EnergyRangeConfig = EnergyRangeConfig()

# class FluxPointsConfig(FeupyBaseConfig):
#     energy: EnergyAxisConfig = EnergyAxisConfig()
#     source: str = "source"
#     parameters: dict = {"selection_optional": "all"}

# class ExcessMapConfig(FeupyBaseConfig):
#     correlation_radius: AngleType = "0.1 deg"
#     parameters: dict = {}
#     energy_edges: EnergyAxisConfig = EnergyAxisConfig()

        
# # Background and Masking Configuration
# class BackgroundConfig(FeupyBaseConfig):
#     method: BackgroundMethodEnum = None
#     exclusion: Path = None
#     parameters: dict = {}


# class SafeMaskConfig(FeupyBaseConfig):
#     methods: List[SafeMaskMethodsEnum] = [SafeMaskMethodsEnum.aeff_default]
#     parameters: dict = {}


# # Dataset and Region Configuration
# class SpatialCircleConfig(FeupyBaseConfig):
#     frame: FrameEnum = None
#     lon: AngleType = None
#     lat: AngleType = None
#     radius: AngleType = None

# class TimeRangeConfig(FeupyBaseConfig):
#     start: TimeType = None
#     stop: TimeType = None
        
# class LightCurveConfig(FeupyBaseConfig):
#     time_intervals: TimeRangeConfig = TimeRangeConfig()
#     energy_edges: EnergyAxisConfig = EnergyAxisConfig()
#     source: str = "source"
#     parameters: dict = {"selection_optional": "all"}
        
        
# class ObservationsConfig(FeupyBaseConfig):
#     datastore: Path = Path("$GAMMAPY_DATA/hess-dl3-dr1/")
#     obs_ids: List[int] = []
#     obs_file: Path = None
#     obs_cone: SpatialCircleConfig = SpatialCircleConfig()
#     obs_time: TimeRangeConfig = TimeRangeConfig()
#     required_irf: List[RequiredHDUEnum] = ["aeff", "edisp", "psf", "bkg"]

# class DatasetsConfig(FeupyBaseConfig):
#     type: ReductionTypeEnum = ReductionTypeEnum.spectrum
#     stack: bool = True
#     geom: GeomConfig = GeomConfig()
#     map_selection: List[MapSelectionEnum] = MapDatasetMaker.available_selection
#     background: BackgroundConfig = BackgroundConfig()
#     safe_mask: SafeMaskConfig = SafeMaskConfig()
#     on_region: SpatialCircleConfig = SpatialCircleConfig()
#     containment_correction: bool = True
#     use_region_center: bool = False
        
# class DatasetOnOffConfig(FeupyBaseConfig):
#     geom: GeomConfig = GeomConfig()
#     map_selection: List[MapSelectionEnum] = MapDatasetMaker.available_selection
#     background: BackgroundConfig = BackgroundConfig()
#     safe_mask: SafeMaskConfig = SafeMaskConfig()
#     on_region: SpatialCircleConfig = SpatialCircleConfig()
#     containment_correction: bool = False
#     containment: float = 0.68
#     use_region_center: bool = False
#     acceptance: int = 1
#     acceptance_off: int = 5

# class StatisticsConfig(FeupyBaseConfig):
#     alpha: float = None
#     n_obs = int = 1
#     wstat: dict = {}
#     fitted_parameters: dict = {}  
        
# class SensitivityConfig(FeupyBaseConfig):
#     gamma_min: int = None
#     n_sigma: int = None
#     bkg_syst_fraction: float = None
#     table_format: TableEnum = 'fits'