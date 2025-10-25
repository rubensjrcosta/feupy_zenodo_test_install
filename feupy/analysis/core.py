# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""ROI and CTAO Analysis classes driving the high level interface API"""

import logging
import html
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion

from gammapy.data import (
    DataStore, Observation, Observations, FixedPointingInfo, PointingMode
)
from gammapy.datasets import (
    Datasets, FluxPointsDataset, MapDataset, SpectrumDataset, SpectrumDatasetOnOff
)
from gammapy.estimators import (
    FluxPoints, SensitivityEstimator, ExcessMapEstimator, FluxPointsEstimator, LightCurveEstimator
)
from gammapy.makers import (
    DatasetsMaker, FoVBackgroundMaker, MapDatasetMaker, ReflectedRegionsBackgroundMaker,
    RingBackgroundMaker, SafeMaskMaker, SpectrumDatasetMaker
)
from gammapy.maps import Map, MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    SkyModel, Models, DatasetModels, FoVBackgroundModel
)
from gammapy.utils.pbar import progress_bar
from gammapy.utils.scripts import make_path

from feupy.analysis.config import ROIAnalysisConfig, CTAOAnalysisConfig
from feupy.analysis.irfs import Irfs
from feupy.catalog.utils import load_catalogs
from feupy.catalog.hawc import get_flux_points_3hwc, get_flux_points_2hwc
from feupy.catalog.fermi import get_flux_points_2PC, get_flux_points_3PC
from feupy.catalog.lhaaso import get_flux_points_1lhaaso
from feupy.catalog.veritas import generate_unique_name
from feupy.sources import Sources, get_catalog_tag
from feupy.utils.coordinates import convert_pos_config_to_skycoord
from feupy.utils.datasets import cut_energy_flux_points_datasets, flux_points_dataset_from_table
from feupy.utils.table import write_tables_csv, write_tables_fits


__all__ = ["ROIAnalysis", "CTAOAnalysis"]

log = logging.getLogger(__name__)

class ROIAnalysis:
    """Config-driven high level simulation interface.

    It is initialized by default with a set of configuration parameters and values declared in
    an internal high level interface model, though the user can also provide configuration
    parameters passed as a nested dictionary at the moment of instantiation. In that case these
    parameters will overwrite the default values of those present in the configuration file.

    Parameters
    ----------
    config : dict or `ROIAnalysisConfig`
        Configuration options following `ROIAnalysisConfig` schema
    """

    def __init__(self, config):
        self.config = config
        self.config.set_logging()
        self.datasets = None
        self.sources = None
        self._has_fp = None
        self.catalog = None

    def _repr_html_(self):
        try:
            return self.to_html()
        except AttributeError:
            return f"<pre>{html.escape(str(self))}</pre>"
            
    @property
    def models(self):
        if not self.datasets:
            raise RuntimeError("No datasets defined. Impossible to set models.")
        return self.datasets.models

    @models.setter
    def models(self, models):
        self.set_models(models, extend=False)
        
    @property
    def config(self):
        """Simulation configuration (`ROIAnalysisConfig`)"""
        return self._config

    @config.setter
    def config(self, value):
        if isinstance(value, dict):
            self._config = ROIAnalysisConfig(**value)
        elif isinstance(value, ROIAnalysisConfig):
            self._config = value
        else:
            raise TypeError("config must be dict or ROIAnalysisConfig.")

    def _get_catalogs(self):
        """Load source catalogs and filter sources within the ROI."""
        catalogs = []
        config_settings = self.config
        all_catalogs = load_catalogs()
        position = convert_pos_config_to_skycoord(config_settings.roi)
        radius = config_settings.roi.radius

        for catalog in all_catalogs:
            separation = position.separation(catalog.positions)
            mask_roi = separation < radius

            if len(catalog[mask_roi].table):
                catalogs.append(catalog[mask_roi])
        return catalogs

    def _get_sources(self):
        """Retrieve sources from catalogs within the ROI."""
        sources = []
        for catalog in self._get_catalogs():
            sources.extend(catalog)
        self.sources =  Sources(sources)

    def _get_catalog_roi(self):
        """
        Compute separations of sources from the target position and return an Astropy Table.

        Returns
        -------
        result_table : `~astropy.table.Table`
            Table with columns for source name, catalog, RA, Dec, and separation from target.
        """
        names, ras, decs, separations, catalogs, has_fp = [], [], [], [], [], []
        sources = self.sources
        config_settings = self.config
        target_position = convert_pos_config_to_skycoord(config_settings.roi)

        for source in sources:
            source_name = source.name
            source_pos = source.position
            ra, dec = source_pos.ra.deg, source_pos.dec.deg
            catalog = get_catalog_tag(source)
            sep = source_pos.separation(target_position).deg

            names.append(source_name)
            ras.append(ra)
            decs.append(dec)
            catalogs.append(catalog)
            separations.append(sep)
            has_fp.append(self._has_fp[source.name])
            
        result_table = Table()
        result_table['index'] =  [_ for _ in range(0, (len(names)))]
        result_table['source_name'] = names
        result_table['source_label'] = sources.labels
        result_table['catalog'] = catalogs
        result_table['ra'] = ras * u.deg
        result_table['dec'] = decs * u.deg
        result_table['separation'] = separations * u.deg
        result_table['has_fp'] = has_fp
        
        for column in result_table.colnames:
            if column.startswith(("ra", "dec", "sep")):
                result_table[column].format = ".3f"
                
        # Add meta (columns descriptions)
        result_table.meta['description'] = {
            'index': "Unique identifier for each source.",
            'source_name': "Name of the source as listed in the catalog.",
            'source_label': "User-defined label for the source, often used for plotting.",
            'has_fp': "Boolean flag indicating if flux points are available (True/False).",
            'catalog': "Catalog from which the source originates.",
            'ra': "Right Ascension (RA) of the source in degrees, formatted to 3 decimal places.",
            'dec': "Declination (Dec) of the source in degrees, formatted to 3 decimal places.",
            'separation': "Angular separation from a reference position in degrees, formatted to 3 decimal places."
        }
        self.catalog = result_table      


    def run(self):
        """Run the ROI analysis, initializing datasets."""
        self._get_sources()
        self._get_flux_points_datasets()
        self._get_catalog_roi()

    def _get_flux_points_datasets(self):
        """Generate FluxPointsDataset objects from sources within the ROI."""
        models = Models()
        datasets = Datasets()
        dict_fp = {}
        sources = self.sources
        log.info("Generating FluxPointsDataset.")

        for index, source in enumerate(sources):
            label = sources.labels[index]
            tag = get_catalog_tag(source)
            if tag == "psrcat":
                dict_fp[source.name] = False
                log.info("Flux points unavailable for ATNF Pulsar Catalog.")
                continue     
                
            try:
                which = None

                if tag == "2PC":
                    spectral_model = self._create_spectral_model(source, which)
                    flux_points = get_flux_points_2PC(source)
                    dict_fp[source.name] = True
                    self._create_flux_points_dataset(f"PSR {label}", models, datasets,  spectral_model, flux_points)

                elif tag == "3PC":
                    spectral_model = self._create_spectral_model(source, which)
                    flux_points = get_flux_points_3PC(source)
                    dict_fp[source.name] = True
                    # self._create_flux_points_dataset(f"PSR {label}", models, datasets,  spectral_model, flux_points)
                    self._create_flux_points_dataset(f"{label}", models, datasets,  spectral_model, flux_points)
                    
                elif tag == "3hwc":
                    spectral_model = self._create_spectral_model(source, which)
                    flux_points = get_flux_points_3hwc(source)
                    dict_fp[source.name] = True
                    # self._create_flux_points_dataset(label, models, datasets,  spectral_model, flux_points)
                    self._create_flux_points_dataset(f"{label}", models, datasets,  spectral_model, flux_points)

                elif tag == "2hwc":
                    for model_type in ['point', 'extended']:
                        spectral_model = self._create_spectral_model(source, model_type)
                        flux_points = get_flux_points_2hwc(source, model_type)
                        dict_fp[source.name] = True
                        self._create_flux_points_dataset(f"{label} ({model_type})", models, datasets,  spectral_model, flux_points)
                         
                elif tag == "1LHAASO":
                    for model_type in ['KM2A', 'WCDA']:
                        spectral_model = self._create_spectral_model(source, model_type)
                        flux_points = get_flux_points_1lhaaso(source, model_type)
                        dict_fp[source.name] = True
                        self._create_flux_points_dataset(f"{label} ({model_type})", models, datasets,  spectral_model, flux_points)
                
                elif tag == "vtscat":
                    unique_names = [source.name]
                    for flux_points in source.flux_points():
                        reference_id = flux_points.meta['reference_id']
                        model_name = generate_unique_name(source.name, reference_id, unique_names)
                        unique_names.append(model_name)
                        model = flux_points.reference_model
                        dict_fp[source.name] = True
                        spectral_model = model.spectral_model
                        label = model_name
                        self._create_flux_points_dataset(label, models, datasets,  spectral_model, flux_points)
                        
                else:        
                    spectral_model = source.spectral_model()
                    flux_points = source.flux_points
                    dict_fp[source.name] = True
                    self._create_flux_points_dataset(label, models, datasets,  spectral_model, flux_points)

                        
            except Exception as error:
                dict_fp[source.name] = False
                log.info(f"Unable to create dataset for {source.name}. Error: {error}")

        self.datasets = datasets
        self.datasets.models = models
        self._has_fp = dict_fp

    @staticmethod
    def _create_spectral_model(source, which=None):
        """Helper function to create a FluxPointsDataset for a source."""
        return source.spectral_model(which=which) if which else source.spectral_model()

    def _create_flux_points_dataset(self,  name, models, datasets,  spectral_model, flux_points):
        """Helper function to create a FluxPointsDataset for a source."""
        
        e_ref_min = self.config.energy_range.min
        e_ref_max = self.config.energy_range.max
        log.info("Getting FluxPointsDataset.")
        
        model = SkyModel(
        name=name,
        spectral_model=spectral_model,
        datasets_names=name
        )
        dataset = FluxPointsDataset(
        models=model,
        data=flux_points,
        name=name
        )
        
        if any([e_ref_min !=  None, e_ref_max !=  None]):
            dataset = cut_energy_flux_points_datasets(
            dataset, 
            e_ref_min, 
            e_ref_max
            ) 
                    
        models.append(model)
        datasets.append(dataset)
        log.info(f"Dataset created: {dataset}")

    def set_models(self, models, extend=True):
        """Set models on datasets.
        Adds `FoVBackgroundModel` if not present already

        Parameters
        ----------
        models : `~gammapy.modeling.models.Models` or str
            Models object or YAML models string
        extend : bool
            Extend the exiting models on the datasets or replace them.
        """
        if not self.datasets or len(self.datasets) == 0:
            raise RuntimeError("Missing datasets")

        log.info("Reading model.")
        if isinstance(models, str):
            models = Models.from_yaml(models)
        elif isinstance(models, Models):
            pass
        elif isinstance(models, DatasetModels) or isinstance(models, list):
            models = Models(models)
        else:
            raise TypeError(f"Invalid type: {models!r}")

        if extend:
            models.extend(self.datasets.models)

        self.datasets.models = models
        
        log.info(models)

    def read_models(self, path, extend=True):
        """Read models from YAML file.

        Parameters
        ----------
        path : str
        Path to the model file.
        extend : bool, optional
        Extend the exiting models on the datasets or replace them.
        Default is True.
        """
        path = make_path(path)
        models = Models.read(path)
        self.set_models(models, extend=extend)
        log.info(f"Models loaded from {path}.")
 
    def read_sources(self, path):
        """Read sources from YAML file.

        Parameters
        ----------
        path : str
        Path to the model file.
        extend : bool, optional
        Extend the exiting sources on the datasets or replace them.
        Default is True.
        """
        path = make_path(path)
        _sources = Sources()
        _sources.read(path)
        self.sources = Sources(_sources)
        log.info(f"sources loaded from {path}.")
    
    def write_sources(self, overwrite=True):
        """Write sources to YAML file.

        File name is taken from the configuration file.
        """
        filename_sources = self.config.general.sources_file
        if filename_sources is not None:
            self.sources.write(
                filename_sources, overwrite=overwrite)
            log.info(f"sources loaded from {filename_sources}.")
        else:
            raise RuntimeError("Missing sources_file in config.general")
            
    def write_models(self, overwrite=True, write_covariance=True):
        """Write models to YAML file.

        File name is taken from the configuration file.
        """
        filename_models = self.config.general.models_file
        if filename_models is not None:
            self.models.write(
                filename_models, overwrite=overwrite, write_covariance=write_covariance
            )
            log.info(f"Models loaded from {filename_models}.")
        else:
            raise RuntimeError("Missing models_file in config.general")

    def read_datasets(self):
        """Read datasets from YAML file.

        File names are taken from the configuration file.
        """
        filename = self.config.general.datasets_file
        filename_models = self.config.general.models_file
        if filename is not None:
            self.datasets = Datasets.read(filename)
            log.info(f"Datasets loaded from {filename}.")
            if filename_models is not None:
                self.read_models(filename_models, extend=False)
        else:
            raise RuntimeError("Missing datasets_file in config.general")

    def write_datasets(self, overwrite=True, write_covariance=True):
        """Write datasets to YAML file.

        File names are taken from the configuration file.

        Parameters
        ----------
        overwrite : bool, optional
            Overwrite existing file. Default is True.
        write_covariance : bool, optional
            Save covariance or not. Default is True.
        """
        filename = self.config.general.datasets_file
        filename_models = self.config.general.models_file
        if filename is not None:
            self.datasets.write(
                filename,
                filename_models,
                overwrite=overwrite,
                write_covariance=write_covariance,
            )
            log.info(f"Datasets stored to {filename}.")
            log.info(f"Datasets stored to {filename_models}.")
        else:
            raise RuntimeError("Missing datasets_file in config.general")


class CTAOAnalysis:
    """Config-driven high level analysis interface.

    It is initialized by default with a set of configuration parameters and values declared in
    an internal high level interface model, though the user can also provide configuration
    parameters passed as a nested dictionary at the moment of instantiation. In that case these
    parameters will overwrite the default values of those present in the configuration file.

    Parameters
    ----------
    config : dict or `~gammapy.analysis.AnalysisConfig`
        Configuration options following `AnalysisConfig` schema.
    """

    def __init__(self, config):
        self.config = config
        self.config.set_logging()
        self.datastore = None
        self.observations = None
        self.datasets = None
        self.fit = Fit()
        self.fit_result = None
        self.flux_points = None

    def _repr_html_(self):
        try:
            return self.to_html()
        except AttributeError:
            return f"<pre>{html.escape(str(self))}</pre>"

    @property
    def models(self):
        if not self.datasets:
            raise RuntimeError("No datasets defined. Impossible to set models.")
        return self.datasets.models

    @models.setter
    def models(self, models):
        self.set_models(models, extend=False)

    @property
    def config(self):
        """Analysis configuration as an `~gammapy.analysis.AnalysisConfig` object."""
        return self._config

    @config.setter
    def config(self, value):
        if isinstance(value, dict):
            self._config = CTAOAnalysisConfig(**value)
        elif isinstance(value, CTAOAnalysisConfig):
            self._config = value
        else:
            raise TypeError("config must be dict or CTAOAnalysisConfig.")

    def simulate_observation(self, obs_id=0):
        """
        Simulate observation with given parameters

        Parameters
        ----------
        obs_settings : `~gammapy.scripts.ObservationParameters`
            Observation parameters
        """
                
        if not self.observations:
            observations = Observations()
        else: 
            observations = self.observations

        if obs_id in observations.ids:
                raise (ValueError("Observation ids must be unique"))
                
        on_region_settings = self.config.datasets.on_region
        observation_settings = self.config.observation
        
        on_lon = on_region_settings.lon
        on_lat = on_region_settings.lat
        frame = on_region_settings.frame
        log.info("Creating the pointing.")
        on_center = SkyCoord(on_lon, on_lat, frame=frame)
        log.info("\n ON center:\n{}".format(on_center))
        
        position_angle = observation_settings.position_angle
        separation = observation_settings.offset
        log.info("\n Obsevation offset:\n{}".format(separation))
        pointing_position = self._create_pointing_position(on_center, position_angle, separation)
        log.info(f"\nPointing position:\n{pointing_position}\n")
        pointing = self._create_pointing(pointing_position)
        log.info(f"\nPointing:\n{pointing}\n")
        
        log.info("\nSetting observation parameters.")
        livetime = observation_settings.livetime
        required_irfs = observation_settings.required_irfs
        irfs = Irfs.get_irfs(required_irfs)
        location = Irfs.get_obs_loc(required_irfs)
        log.info("\nirfs: {}".format(Irfs.get_irfs_label(required_irfs)))
        log.info("\nlocation: {}".format(location))
        log.info("\nlivetime: {}".format(livetime))
        observation = Observation.create(pointing=pointing, livetime=livetime, irfs=irfs, location=location, obs_id=obs_id)
                        
        observations.append(observation)
        self.observations = observations
        log.info(f"\n{observation}\n")
        log.info(f"Observation {observation.obs_id} loaded.")

    
    @staticmethod
    def _create_pointing_position(position, position_angle, separation):
        """Create the pointing position"""
        return position.directional_offset_by(position_angle, separation)

    @staticmethod
    def _create_pointing(pointing_position):
        """Create the pointing."""
        return FixedPointingInfo(
            mode=PointingMode.POINTING,
            fixed_icrs=pointing_position.icrs,
    
        )
        
    def get_spectrum_dataset(self, model=None, obs_id=0, random_state=42):
        """
        Produce reduced datasets.

        Notes
        -----
        The progress bar can be displayed for this function.
        """
        datasets_settings = self.config.datasets
        if not self.observations or len(self.observations) == 0:
            raise RuntimeError("No observations have been selected.")

        if datasets_settings.type == "1d":
            self._spectrum_extraction(model=model, obs_id=obs_id, random_state=random_state)
        else: raise ValueError(
                    f"Incorrect dataset type. Expect '1d'. Got {datasets_settings.type}."
                )
    def update_config(self, config):
        """Update the configuration."""
        self.config = self.config.update(config=config)
        
    def _spectrum_extraction(self, model=None, obs_id=0, random_state=42):
        """Make the SpectrumDataset for ON-OFF analysis"""
        datasets_settings = self.config.datasets
        obs_settings = self.config.observation   
        energy_axis = self._make_energy_axis(self.config.datasets.geom.axes.energy)
        energy_axis_true = self._make_energy_axis(self.config.datasets.geom.axes.energy_true)

        dataset_maker = self._create_dataset_maker()
        safe_mask_maker = self._create_safe_mask_maker()
        bkg_maker = self._create_background_maker()
        
        log.info("Getting the observation.")
        observation = self.observations[obs_id]
        log.info(f"\n{observation}\n")
        log.info("Getting the reference Dataset.")    
        reference = self._create_reference_dataset(str(observation.obs_id))
        log.info("\nreference: {}".format(reference))
        
        log.info("Reducing spectrum datasets.")
        dataset = dataset_maker.run(reference, observation)
        log.info("\nMaker: {}".format(dataset))

        if not datasets_settings.containment_correction:
            
            if datasets_settings.containment:
                
                containment = datasets_settings.containment

                dataset.exposure *= containment
                log.info("\nCorrected exposure (containment: {}%):\n{}\n".format(containment, dataset)) 

                offset = obs_settings.offset

                on_radii = observation.psf.containment_radius(
                    energy_true=energy_axis.center, offset=offset, fraction=containment
                )
                self._on_radii = on_radii
                on_region_radius = datasets_settings.on_region.radius
                factor = (1 - np.cos(on_radii)) / (1 - np.cos(on_region_radius))
                dataset.background *= factor.value.reshape((-1, 1, 1))
                log.info("\nCorrected background (containment: {}%):\n{}\n".format(containment, dataset)) 
                
        if bkg_maker is not None:
            dataset = bkg_maker.run(dataset, observation)
            if dataset.counts_off is None:
                raise ValueError(
                    f"No OFF region found for observation {observation.obs_id}. Discarding."
                )
        dataset = safe_mask_maker.run(dataset, observation)
        log.info("\nSafe mask maker: {}".format(dataset))
        
        
        if model is not None:
            dataset.models = model
            dataset.fake(random_state=random_state)
            log.info(f"\nDataset Model:\n{dataset}\n")      
        
        self.spectrum_dataset = dataset
        log.info(f"\nDataset:\n{self.spectrum_dataset}\n")

    def _create_dataset_maker(self):
        """Create the Dataset Maker."""
        log.debug("Creating the target Dataset Maker.")

        datasets_settings = self.config.datasets
        if datasets_settings.type == "3d":
            maker = MapDatasetMaker(selection=datasets_settings.map_selection)
        elif datasets_settings.type == "1d":
            maker_config = {}
            if datasets_settings.containment_correction:
                maker_config["containment_correction"] = (
                    datasets_settings.containment_correction
                )

            maker_config["selection"] = datasets_settings.map_selection
            maker_config["use_region_center"] = datasets_settings.use_region_center
            maker = SpectrumDatasetMaker(**maker_config)

        return maker

    def _create_safe_mask_maker(self):
        """Create the SafeMaskMaker."""
        log.debug("Creating the mask_safe Maker.")

        safe_mask_selection = self.config.datasets.safe_mask.methods
        safe_mask_settings = self.config.datasets.safe_mask.parameters
        return SafeMaskMaker(methods=safe_mask_selection, **safe_mask_settings)

    def _create_background_maker(self):
        """Create the Background maker."""
        log.info("Creating the background Maker.")

        datasets_settings = self.config.datasets
        bkg_maker_config = {}
        if datasets_settings.background.exclusion:
            path = make_path(datasets_settings.background.exclusion)
            exclusion_mask = Map.read(path)
            exclusion_mask.data = exclusion_mask.data.astype(bool)
            bkg_maker_config["exclusion_mask"] = exclusion_mask
        bkg_maker_config.update(datasets_settings.background.parameters)

        bkg_method = datasets_settings.background.method

        bkg_maker = None
        if bkg_method == "fov_background":
            log.debug(f"Creating FoVBackgroundMaker with arguments {bkg_maker_config}")
            bkg_maker = FoVBackgroundMaker(**bkg_maker_config)
        elif bkg_method == "ring":
            bkg_maker = RingBackgroundMaker(**bkg_maker_config)
            log.debug(f"Creating RingBackgroundMaker with arguments {bkg_maker_config}")
            if datasets_settings.geom.axes.energy.nbins > 1:
                raise ValueError(
                    "You need to define a single-bin energy geometry for your dataset."
                )
        elif bkg_method == "reflected":
            bkg_maker = ReflectedRegionsBackgroundMaker(**bkg_maker_config)
            log.debug(
                f"Creating ReflectedRegionsBackgroundMaker with arguments {bkg_maker_config}"
            )
        else:
            log.warning("No background maker set. Check configuration.")
        return bkg_maker

    def _create_reference_dataset(self, name=None):
        """Create the reference dataset for the current analysis."""
        log.debug("Creating target Dataset.")
        geom = self._create_geometry()

        geom_settings = self.config.datasets.geom
        geom_irf = dict(energy_axis_true=None, binsz_irf=None)
        if geom_settings.axes.energy_true.min is not None:
            geom_irf["energy_axis_true"] = self._make_energy_axis(
                geom_settings.axes.energy_true, name="energy_true"
            )
        if geom_settings.wcs.binsize_irf is not None:
            geom_irf["binsz_irf"] = geom_settings.wcs.binsize_irf.to("deg").value

        if self.config.datasets.type == "1d":
            return SpectrumDataset.create(geom, name=name, **geom_irf)
        else:
            return MapDataset.create(geom, name=name, **geom_irf)

    @staticmethod
    def _create_region_geometry(on_region_settings, axes):
        """Create the region geometry."""
        on_lon = on_region_settings.lon
        on_lat = on_region_settings.lat
        on_center = SkyCoord(on_lon, on_lat, frame=on_region_settings.frame)
        on_region = CircleSkyRegion(on_center, on_region_settings.radius)

        return RegionGeom.create(region=on_region, axes=axes)

    def _create_geometry(self):
        """Create the geometry."""
        log.debug("Creating geometry.")
        datasets_settings = self.config.datasets
        geom_settings = datasets_settings.geom
        axes = [self._make_energy_axis(geom_settings.axes.energy)]
        if datasets_settings.type == "3d":
            geom = self._create_wcs_geometry(geom_settings.wcs, axes)
        elif datasets_settings.type == "1d":
            geom = self._create_region_geometry(datasets_settings.on_region, axes)
        else:
            raise ValueError(
                f"Incorrect dataset type. Expect '1d' or '3d'. Got {datasets_settings.type}."
            )
        return geom

    @staticmethod
    def _create_wcs_geometry(wcs_geom_settings, axes):
        """Create the WCS geometry."""
        geom_params = {}
        skydir_settings = wcs_geom_settings.skydir
        if skydir_settings.lon is not None:
            skydir = SkyCoord(
                skydir_settings.lon, skydir_settings.lat, frame=skydir_settings.frame
            )
            geom_params["skydir"] = skydir

        if skydir_settings.frame in ["icrs", "galactic"]:
            geom_params["frame"] = skydir_settings.frame
        else:
            raise ValueError(
                f"Incorrect skydir frame: expect 'icrs' or 'galactic'. Got {skydir_settings.frame}"
            )

        geom_params["axes"] = axes
        geom_params["binsz"] = wcs_geom_settings.binsize
        width = wcs_geom_settings.width.width.to("deg").value
        height = wcs_geom_settings.width.height.to("deg").value
        geom_params["width"] = (width, height)

        return WcsGeom.create(**geom_params)
        
    def get_datasets(self):
        """
        Produce reduced datasets.

        Notes
        -----
        The progress bar can be displayed for this function.
        """
        datasets_settings = self.config.datasets
        
        if not self.observations or len(self.observations) == 0:
            raise RuntimeError("No observations have been selected.")

        if not self.spectrum_dataset:
            raise RuntimeError("No spectrum dataset have been selected.")

        
        if datasets_settings.type == "1d":
            self._run_on_off()
            
        else: raise ValueError(
                    f"Incorrect dataset type. Expect '1d'. Got {datasets_settings.type}."
                )

    def run_fit(self):
        """Fitting reduced datasets to model."""
        if not self.models:
            raise RuntimeError("Missing models")

        fit_settings = self.config.fit
        for dataset in self.datasets:
            if fit_settings.fit_range:
                energy_min = fit_settings.fit_range.min
                energy_max = fit_settings.fit_range.max
                geom = dataset.counts.geom
                dataset.mask_fit = geom.energy_mask(energy_min, energy_max)

        log.info("Fitting datasets.")
        result = self.fit.run(datasets=self.datasets)
        self.fit_result = result
        log.info(self.fit_result)
        
    def _run_on_off(self):
        datasets_settings = self.config.datasets
        # on_off_settings = self.config.on_off
        stat_settings = self.config.statistics
        
        dataset = self.spectrum_dataset
        
        acceptance = datasets_settings.on_off.acceptance 
        acceptance_off = datasets_settings.on_off.acceptance_off
        
        dataset_on_off = self._create_dataset_on_off(dataset, acceptance, acceptance_off)
        
        n_obs =  stat_settings.n_obs
        datasets = Datasets()
        for idx in range(n_obs):
            dataset_on_off.fake(random_state=idx, npred_background=dataset.npred_background())
            dataset_fake = dataset_on_off.copy(name=f"obs-{idx}")
            dataset_fake.meta_table["OBS_ID"] = [idx]
            datasets.append(dataset_fake)
        table = datasets.info_table()
        display(table)
        # self.datasets_on_off = datasets
        self.datasets = datasets

        if datasets_settings.stack:
            stacked = self.datasets.stack_reduce(name="stacked")
            self.datasets = Datasets([stacked])
    
    def get_flux_points(self):
        """Calculate flux points for a specific model component."""
        if not self.datasets:
            raise RuntimeError(
                "No datasets defined. Impossible to compute flux points."
            )

        fp_settings = self.config.flux_points
        log.info("Calculating flux points.")
        energy_edges = self._make_energy_axis(fp_settings.energy).edges
        flux_point_estimator = FluxPointsEstimator(
            energy_edges=energy_edges,
            source=fp_settings.source,
            fit=self.fit,
            n_jobs=self.config.general.n_jobs,
            **fp_settings.parameters,
        )

        fp = flux_point_estimator.run(datasets=self.datasets)

        self.flux_points = FluxPointsDataset(
            data=fp, models=self.models[fp_settings.source]
        )
        cols = ["e_ref", "dnde", "dnde_ul", "dnde_err", "sqrt_ts"]
        table = self.flux_points.data.to_table(sed_type="dnde")
        log.info("\n{}".format(table[cols]))                

    @staticmethod    
    def _create_dataset_on_off(dataset, acceptance, acceptance_off):
    # Spectrum dataset for on-off likelihood fitting.
        dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
            dataset=dataset, 
            acceptance=acceptance, 
            acceptance_off=acceptance_off,
        )
        dataset_on_off.fake(
            random_state='random-seed', 
            npred_background=dataset.npred_background()
        )
        return(dataset_on_off)
        
    @staticmethod
    def _make_energy_axis(axis, name="energy"):
        if axis.min is None or axis.max is None:
            return None
        elif axis.nbins is None or axis.nbins < 1:
            return None
        else:
            return MapAxis.from_bounds(
                name=name,
                lo_bnd=axis.min.value,
                hi_bnd=axis.max.to_value(axis.min.unit),
                nbin=axis.nbins,
                unit=axis.min.unit,
                interp="log",
                node_type="edges",
            )

    def set_models(self, models, extend=True):
        """Set models on datasets.

        Adds `FoVBackgroundModel` if not present already

        Parameters
        ----------
        models : `~gammapy.modeling.models.Models` or str
            Models object or YAML models string.
        extend : bool, optional
            Extend the exiting models on the datasets or replace them.
            Default is True.
        """
        if not self.datasets or len(self.datasets) == 0:
            raise RuntimeError("Missing datasets")

        log.info("Reading model.")
        if isinstance(models, str):
            models = Models.from_yaml(models)
        elif isinstance(models, Models):
            pass
        elif isinstance(models, DatasetModels) or isinstance(models, list):
            models = Models(models)
        else:
            raise TypeError(f"Invalid type: {models!r}")

        if extend:
            models.extend(self.datasets.models)

        self.datasets.models = models

        bkg_models = []
        for dataset in self.datasets:
            if dataset.tag == "MapDataset" and dataset.background_model is None:
                bkg_models.append(FoVBackgroundModel(dataset_name=dataset.name))
        if bkg_models:
            models.extend(bkg_models)
            self.datasets.models = models

        log.info(models)

    def read_models(self, path, extend=True):
        """Read models from YAML file.

        Parameters
        ----------
        path : str
            Path to the model file.
        extend : bool, optional
            Extend the exiting models on the datasets or replace them.
            Default is True.
        """
        path = make_path(path)
        models = Models.read(path)
        self.set_models(models, extend=extend)
        log.info(f"Models loaded from {path}.")

    def write_models(self, overwrite=True, write_covariance=True):
        """Write models to YAML file.

        File name is taken from the configuration file.
        """
        filename_models = self.config.general.models_file
        if filename_models is not None:
            self.models.write(
                filename_models, overwrite=overwrite, write_covariance=write_covariance
            )
            log.info(f"Models loaded from {filename_models}.")
        else:
            raise RuntimeError("Missing models_file in config.general")

    def read_datasets(self):
        """Read datasets from YAML file.

        File names are taken from the configuration file.
        """
        filename = self.config.general.datasets_file
        filename_models = self.config.general.models_file
        if filename is not None:
            self.datasets = Datasets.read(filename)
            log.info(f"Datasets loaded from {filename}.")
            if filename_models is not None:
                self.read_models(filename_models, extend=False)
        else:
            raise RuntimeError("Missing datasets_file in config.general")

    def write_datasets(self, overwrite=True, write_covariance=True):
        """Write datasets to YAML file.

        File names are taken from the configuration file.

        Parameters
        ----------
        overwrite : bool, optional
            Overwrite existing file. Default is True.
        write_covariance : bool, optional
            Save covariance or not. Default is True.
        """
        filename = self.config.general.datasets_file
        filename_models = self.config.general.models_file
        if filename is not None:
            self.datasets.write(
                filename,
                filename_models,
                overwrite=overwrite,
                write_covariance=write_covariance,
            )
            log.info(f"Datasets stored to {filename}.")
            log.info(f"Datasets stored to {filename_models}.")
        else:
            raise RuntimeError("Missing datasets_file in config.general")


    # Sensitivity Computation
    def compute_sensitivity(self):
        """Compute the sensitivity of the analysis."""
        sens_settings = self.config.sensitivity
        obs_settings = self.config.observation
        energy_settings = self.config.datasets.geom.axes.energy
        datasets_settings = self.config.datasets
        dataset = self.spectrum_dataset

        acceptance = datasets_settings.on_off.acceptance 
        acceptance_off = datasets_settings.on_off.acceptance_off
        
        


        sensitivity_estimator = SensitivityEstimator(
            gamma_min=sens_settings.gamma_min,
            n_sigma=sens_settings.n_sigma,
            bkg_syst_fraction=sens_settings.bkg_syst_fraction
        )
        
        log.info("\nCreating ON/OFF Dataset:")
        dataset_on_off = self._create_dataset_on_off(dataset, acceptance, acceptance_off)
        log.info("\nON/OFF Dataset:\n{}\n".format(dataset_on_off))
        
        table = sensitivity_estimator.run(dataset_on_off)

        if self._on_radii is not None:
            table["on_radii"] = self._on_radii
            table["on_radii"].format = '.3e'

        required_irfs = obs_settings.required_irfs
        irfs_label = Irfs.get_irfs_label(required_irfs, which='both')
        dataset_label = f'sens {irfs_label}'

        dataset = flux_points_dataset_from_table(table, name=dataset_label)
        self.dataset_sens = dataset

        dataset_on_off1 = dataset_on_off.to_image()

        sensitivity_estimator1 = SensitivityEstimator(
            gamma_min=sens_settings.gamma_min,
            n_sigma=sens_settings.n_sigma,
            bkg_syst_fraction=sens_settings.bkg_syst_fraction
        )
        table1 = sensitivity_estimator1.run(dataset_on_off1)
        log.info("\n{}".format(table1))

        # Convert to a `FluxPoints` object for integral flux
        flux_points = FluxPoints.from_table(
            table1,
            sed_type="e2dnde",
            reference_model=sensitivity_estimator1.spectral_model
        )
        int_sens = np.squeeze(flux_points.flux.quantity)

        livetime = obs_settings.livetime
        log.info(
            f"Integral sensitivity in {livetime:.2f} above {energy_settings.min:.2e} "
            f"is {int_sens:.2e}"
        )
        
        table.meta = self._get_table_meta()
        table.meta['INT_SENS'] = f"{int_sens:.2e}"
        self.table_sens = table

    def _get_table_meta(self):
        """Retrieve metadata for the sensitivity table."""
        obs_settings = self.config.observation
        meta = {
            # "SOURCE": obs_settings.target.name,
            "ONRADIUS": f"{self.config.datasets.on_region.radius.to('deg').value} deg",
            "OFFSET": obs_settings.offset.to_string(),
            "LIVETIME": obs_settings.livetime.to_string(),
            "ARRAY": Irfs.get_irfs_array(obs_settings.required_irfs),
            "AZIMUTH": obs_settings.required_irfs[1],
            "ZENITH": u.Quantity(obs_settings.required_irfs[2]).to_string(),
            "OBS_TIME": u.Quantity(obs_settings.required_irfs[3]).to_string(),
            'IRFS': obs_settings.required_irfs
        }
        return meta

    # File Handling Methods
    def write_table_sensitivity(self, overwrite=True):
        """Write sensitivity table to file in specified format."""
        if self.table_sens is not None:
            path_file = self.config.sensitivity.data_path
            file_name = self.get_file_name()
            file_format = self.config.sensitivity.table_format
            
            if file_format == 'csv':
                write_tables_csv(self.table_sens, path_file, file_name)
                log.info(f"Table ({file_name}.csv) stored to {path_file}.")
            else:
                write_tables_fits(self.table_sens, path_file, file_name)
                log.info(f"Table ({file_name}.fits) stored to {path_file}.")
        else:
            raise RuntimeError("Missing table_sens")

    def read_table_sensitivity(self):
        """Read datasets from file based on the configuration format."""
        try:
            path_file = self.config.sensitivity.data_path
            file_name = self.get_file_name()
            if self.config.sensitivity.table_format == 'csv':
                return Table.read(f'{path_file}/{file_name}.csv', format='ascii')
            else:
                return Table.read(f'{path_file}/{file_name}.fits', format='fits')
        except Exception as error:
            log.error(f'Error reading sensitivity table: {error}')
            raise

    def get_file_name(self, which='both'):
        """Generate file name for saving sensitivity data."""
        obs_settings = self.config.observation
        irfs_name = Irfs.get_irfs_name(obs_settings.required_irfs, which=which)
        dataset_name = f'sens_{irfs_name}'
        livetime = obs_settings.livetime        
        return f"{dataset_name.replace(' ', '-')}_livetime{livetime.to_string().replace(' ', '')}"


