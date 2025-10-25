# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""ROI and CTAO Analysis configuration."""

import json
import logging
from collections import defaultdict
from pathlib import Path

import yaml

from gammapy.utils.scripts import make_path, read_yaml
from gammapy.analysis.config import (
    GammapyBaseConfig, GeneralConfig, EnergyRangeConfig, SpatialCircleConfig,
    deep_update, FluxPointsConfig, FitConfig, SkyCoordConfig
)

from feupy.utils.config import (
    ObservationConfig, DatasetsConfig, OnOffConfig, SensitivityConfig, StatisticsConfig
)


__all__ = ["ROIAnalysisConfig", "CTAOAnalysisConfig"]

CONFIG_PATH = Path(__file__).resolve().parent / "config"
DOCS_FILE = CONFIG_PATH / "docs.yaml"

log = logging.getLogger(__name__)

        
class ROIAnalysisConfig(GammapyBaseConfig):
    """Gammapy analysis configuration."""

    general: GeneralConfig = GeneralConfig()
    roi: SpatialCircleConfig = SpatialCircleConfig()
    energy_range: EnergyRangeConfig = EnergyRangeConfig()

    def __str__(self):
        """Display settings in pretty YAML format."""
        info = self.__class__.__name__ + "\n\n\t"
        data = self.to_yaml()
        data = data.replace("\n", "\n\t")
        info += data
        return info.expandtabs(tabsize=4)

    @classmethod
    def read(cls, path):
        """Read from YAML file."""
        config = read_yaml(path)
        return ROIAnalysisConfig(**config)

    @classmethod
    def from_yaml(cls, config_str):
        """Create from YAML string."""
        settings = yaml.safe_load(config_str)
        return ROIAnalysisConfig(**settings)

    def write(self, path, overwrite=False):
        """Write to YAML file."""
        path = make_path(path)

        if path.exists() and not overwrite:
            raise IOError(f"File exists already: {path}")

        path.write_text(self.to_yaml())

    def to_yaml(self):
        """Convert to YAML string."""
        data = json.loads(self.model_dump_json())
        return yaml.dump(
            data, sort_keys=False, indent=4, width=80, default_flow_style=None
        )

    def set_logging(self):
        """Set logging config.

        Calls ``logging.basicConfig``, i.e. adjusts global logging state.
        """
        self.general.log.level = self.general.log.level.upper()
        logging.basicConfig(**self.general.log.model_dump())
        log.info("Setting logging config: {!r}".format(self.general.log.model_dump()))

    def update(self, config=None):
        """Update config with provided settings.

        Parameters
        ----------
        config : str or `ROIAnalysisConfig` object, optional
            Configuration settings provided in dict() syntax. Default is None.
        """
        if isinstance(config, str):
            other = ROIAnalysisConfig.from_yaml(config)
        elif isinstance(config, ROIAnalysisConfig):
            other = config
        else:
            raise TypeError(f"Invalid type: {config}")

        config_new = deep_update(
            self.model_dump(exclude_defaults=True),
            other.model_dump(exclude_defaults=True),
        )
        return ROIAnalysisConfig(**config_new)

    @staticmethod
    def _get_doc_sections():
        """Return dictionary with commented docs from docs file."""
        doc = defaultdict(str)
        with open(DOCS_FILE) as f:
            for line in filter(lambda line: not line.startswith("---"), f):
                line = line.strip("\n")
                if line.startswith("# Section: "):
                    keyword = line.replace("# Section: ", "")
                doc[keyword] += line + "\n"
        return doc


class CTAOAnalysisConfig(GammapyBaseConfig):
    """Gammapy analysis configuration."""

    general: GeneralConfig = GeneralConfig()
    observation: ObservationConfig = ObservationConfig()
    datasets: DatasetsConfig = DatasetsConfig()
    statistics: StatisticsConfig = StatisticsConfig()
    fit: FitConfig = FitConfig()
    flux_points: FluxPointsConfig = FluxPointsConfig()
    sensitivity: SensitivityConfig = SensitivityConfig()

    # energy_range: EnergyRangeConfig = EnergyRangeConfig()

    def __str__(self):
        """Display settings in pretty YAML format."""
        info = self.__class__.__name__ + "\n\n\t"
        data = self.to_yaml()
        data = data.replace("\n", "\n\t")
        info += data
        return info.expandtabs(tabsize=4)

    @classmethod
    def read(cls, path):
        """Read from YAML file."""
        config = read_yaml(path)
        return CTAOAnalysisConfig(**config)

    @classmethod
    def from_yaml(cls, config_str):
        """Create from YAML string."""
        settings = yaml.safe_load(config_str)
        return CTAOAnalysisConfig(**settings)

    def write(self, path, overwrite=False):
        """Write to YAML file."""
        path = make_path(path)

        if path.exists() and not overwrite:
            raise IOError(f"File exists already: {path}")

        path.write_text(self.to_yaml())

    def to_yaml(self):
        """Convert to YAML string."""
        data = json.loads(self.model_dump_json())
        return yaml.dump(
            data, sort_keys=False, indent=4, width=80, default_flow_style=None
        )

    def set_logging(self):
        """Set logging config.

        Calls ``logging.basicConfig``, i.e. adjusts global logging state.
        """
        self.general.log.level = self.general.log.level.upper()
        logging.basicConfig(**self.general.log.model_dump())
        log.info("Setting logging config: {!r}".format(self.general.log.model_dump()))

    def update(self, config=None):
        """Update config with provided settings.

        Parameters
        ----------
        config : str or `CTAOAnalysisConfig` object, optional
            Configuration settings provided in dict() syntax. Default is None.
        """
        if isinstance(config, str):
            other = CTAOAnalysisConfig.from_yaml(config)
        elif isinstance(config, CTAOAnalysisConfig):
            other = config
        else:
            raise TypeError(f"Invalid type: {config}")

        config_new = deep_update(
            self.model_dump(exclude_defaults=True),
            other.model_dump(exclude_defaults=True),
        )
        return CTAOAnalysisConfig(**config_new)

    @staticmethod
    def _get_doc_sections():
        """Return dictionary with commented docs from docs file."""
        doc = defaultdict(str)
        with open(DOCS_FILE) as f:
            for line in filter(lambda line: not line.startswith("---"), f):
                line = line.strip("\n")
                if line.startswith("# Section: "):
                    keyword = line.replace("# Section: ", "")
                doc[keyword] += line + "\n"
        return doc