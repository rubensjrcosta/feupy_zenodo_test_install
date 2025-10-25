# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Datasets module.

This module provides utilities for handling datasets and flux points datasets, 
including functions for creating, cutting, and retrieving energy bounds for datasets.

Functions:
    - get_energy_bounds: Calculate minimum and maximum energy bounds.
    - cut_energy_flux_points_dataset: Cut flux points dataset within specified energy limits.
    - flux_points_dataset_from_table: Create a `FluxPointsDataset` from a table.
"""

from gammapy.datasets import Datasets, FluxPointsDataset
from gammapy.estimators import FluxPoints
from gammapy.modeling.models import SkyModel
from gammapy.modeling import Fit

from astropy import units as u

import numpy as np

from feupy.utils.scripts import is_documented_by
from feupy.sources import Sources


__all__ = [
    "get_datasets",
    "get_sources",
    "get_energy_bounds_from_datasets",
    "cut_energy_flux_points_datasets",
    "flux_points_dataset_from_table",
]


def get_datasets(datasets, names):
    """
    Retrieve datasets by name from a Datasets object.
    
    Parameters:
    datasets: Datasets
        Collection of datasets.
    names: list of str
        Names of the datasets to retrieve.
        
    Returns:
    Datasets
        Filtered datasets that match the provided names.
    """
    return Datasets([x for x in datasets if x.name in names])


def get_sources(sources, names):
    """
    Retrieve sources by name from a Sources object.
    
    Parameters:
    sources: Sources
        Collection of sources.
    names: list of str
        Names of the sources to retrieve.
        
    Returns:
    Sources
        Filtered sources that match the provided names.
    """
    return Sources([x for x in sources if x.name in names])

def get_energy_bounds_from_datasets(datasets):
    """
    Calculate the energy bounds (minimum and maximum energies) from one or multiple datasets.
    
    Parameters
    ----------
    datasets : `~gammapy.datasets.FluxPointsDataset` or `~gammapy.datasets.Datasets`
        The dataset or collection of datasets containing spectral data with energy bounds.
        
    Returns
    -------
    energy_bounds : `~astropy.units.Quantity`
        A two-element array containing the minimum and maximum energies as astropy quantities.
        If the minimum and maximum energies are equal, a small range is created around the value.
    """
    # Ensure datasets is iterable (handle single dataset case)
    if isinstance(datasets, FluxPointsDataset):
        datasets = [datasets]
    
    # Flatten all energy_min and energy_max values to handle cases where they are arrays
    energy_mins = u.Quantity([min(dataset.data.energy_min).to(u.TeV) for dataset in datasets])
    energy_maxs = u.Quantity([max(dataset.data.energy_max).to(u.TeV) for dataset in datasets])

    energy_min = energy_mins.min()  # Get the minimum energy across all datasets
    energy_max = energy_maxs.max()  # Get the maximum energy across all datasets
    
    if energy_max == energy_min:
        print(f'\nThe minimum and maximum energies are equal, so a small range will be created around the value: {energy_max}')
        energy_max += 1 * energy_max.unit
        energy_min = energy_max - 1 * energy_max.unit
    
    return u.Quantity([energy_min, energy_max])


# def get_energy_bounds(dataset):
#     """
#     Calculate the energy bounds (minimum and maximum energies) from the dataset.
    
#     Parameters
#     ----------
#     dataset : `~gammapy.datasets.Dataset`
#         A dataset containing spectral data with energy bounds.
        
#     Returns
#     -------
#     energy_bounds : list of `~astropy.units.Quantity`
#         A list containing the minimum and maximum energies as astropy quantities.
#         If the minimum and maximum energies are equal, a small range is created around the value.
#     """
#     energy_min = np.min(dataset.data.energy_min)
#     energy_max = np.max(dataset.data.energy_max)
    
#     if energy_max == energy_min:
#         print('\nThe minimum and maximum energies are equal, so a small range will be created around the value: {}'.format(energy_max))
#         energy_max += 1 * energy_max.unit
#         energy_min = energy_max - 1 * energy_max.unit
    
#     return u.Quantity([energy_min, energy_max])

def cut_energy_flux_points_datasets(datasets, e_ref_min=None, e_ref_max=None):
    """
    Cut flux points dataset(s) within specified energy limits.
    
    Parameters
    ----------
    datasets : `~gammapy.datasets.FluxPointsDataset` or `~gammapy.datasets.Datasets`
        The dataset or collection of datasets to be filtered by energy.
    e_ref_min : `~astropy.units.Quantity`, optional
        Minimum energy reference to filter the dataset (inclusive).
    e_ref_max : `~astropy.units.Quantity`, optional
        Maximum energy reference to filter the dataset (inclusive).
        
    Returns
    -------
    `FluxPointsDataset` or `Datasets`
        Filtered dataset(s) with flux points within the specified energy range.
    """
    # If input is a single FluxPointsDataset, wrap it in a list for consistent processing
    single_dataset = False
    if isinstance(datasets, FluxPointsDataset):
        datasets = [datasets]
        single_dataset = True
    
    # Create a new Datasets collection for storing results
    filtered_datasets = Datasets()
    for dataset in datasets:
        try:
            # Skip datasets that are not FluxPointsDataset instances
            if not isinstance(dataset, FluxPointsDataset):
                print(f"Skipping dataset '{dataset.name}' - not a FluxPointsDataset.")
                continue

            flux_points = dataset.data
            models = dataset.models[0] if dataset.models else None
            ds_name = dataset.name
            # Apply minimum energy cut if specified
            if e_ref_min is not None:
                mask_energy = np.array([e_ref >= e_ref_min for e_ref in flux_points.energy_ref])
                flux_points = FluxPoints.from_table(flux_points.to_table()[mask_energy])

            # Apply maximum energy cut if specified
            if e_ref_max is not None:
                mask_energy = np.array([e_ref <= e_ref_max for e_ref in flux_points.energy_ref])
                flux_points = FluxPoints.from_table(flux_points.to_table()[mask_energy])

            # Create and append the filtered FluxPointsDataset
            filtered_dataset = FluxPointsDataset(models=models, data=flux_points, name=ds_name)
            filtered_datasets.append(filtered_dataset)
        except Exception as error:
            print('\nUnable to cut {} FluxPointsDataset. An error has occurred: {}'.format(dataset.name, error))
        # Return a single FluxPointsDataset if input was a single dataset
#     if not filtered_datasets:    
#         filtered_datasets = datasets
    if single_dataset:
        return filtered_datasets[0]
    return filtered_datasets

# def cut_energy_flux_points_dataset(dataset, e_ref_min=None, e_ref_max=None):
#     """
#     Cut flux points dataset within specified energy limits.
    
#     Parameters
#     ----------
#     dataset : `~gammapy.datasets.FluxPointsDataset`
#         The dataset to be filtered by energy.
#     e_ref_min : `~astropy.units.Quantity`, optional
#         Minimum energy reference to filter the dataset (inclusive).
#     e_ref_max : `~astropy.units.Quantity`, optional
#         Maximum energy reference to filter the dataset (inclusive).
        
#     Returns
#     -------
#     FluxPointsDataset
#         Filtered dataset with flux points within the specified energy range.
#     """
#     flux_points = dataset.data
#     models = dataset.models[0] if dataset.models else None
#     ds_name = dataset.name
#     try:
#         if e_ref_min is not None:
#             mask_energy = np.array([e_ref >= e_ref_min for e_ref in flux_points.energy_ref])
#             flux_points = FluxPoints.from_table(flux_points.to_table()[mask_energy])

#         if e_ref_max is not None:
#             mask_energy = np.array([e_ref <= e_ref_max for e_ref in flux_points.energy_ref])
#             flux_points = FluxPoints.from_table(flux_points.to_table()[mask_energy])
#     except Exception as error:
#         print('\nUnable to cut {} FluxPointsDataset. An error has occurred: {}'.format(dataset.name, error))
            
#     return FluxPointsDataset(models=models, data=flux_points, name=ds_name)


@is_documented_by([FluxPoints, FluxPointsDataset])
def flux_points_dataset_from_table(
    table,
    reference_model=None,
    sed_type=None,
    name=None,
    kwargs_fp={'format': 'gadf-sed', 'gti': None},
    kwargs_ds={'mask_fit': None, 'mask_safe': None, 'meta_table': None},
    model_name=None,
):
    """
    Create a `FluxPointsDataset` from a table.
    
    Parameters
    ----------
    table : `~astropy.table.Table`
        Table containing the flux points.
    reference_model : `~gammapy.modeling.models.SpectralModel`, optional
        Reference spectral model.
    sed_type : str, optional
        Type of spectral energy distribution.
    name : str, optional
        Name of the dataset.
    kwargs_fp : dict, optional
        Additional arguments for `FluxPoints.from_table`.
    kwargs_ds : dict, optional
        Additional arguments for `FluxPointsDataset`.
    model_name : str, optional
        Name of the model to use if `reference_model` is provided.
        
    Returns
    -------
    FluxPointsDataset
        Dataset created from the provided table.
    """
    flux_points = FluxPoints.from_table(
        table=table,
        reference_model=reference_model,
        sed_type=sed_type,
        **kwargs_fp,
    )
    
    models = None
    if reference_model:
        models = SkyModel(spectral_model=reference_model, name=model_name or name)
    
    return FluxPointsDataset(
        models=models,
        data=flux_points,
        name=name,
        **kwargs_ds,
    )


# __all__ = [
#     "get_energy_bounds",
#     "flux_points_dataset_from_table",
#     "cut_energy_flux_points_dataset",
#     "cut_energy_flux_points_datasets",
#     "create_spectrum_dataset_empty",
#     "create_spectrum_dataset_onoff"
# ]


# class DatasetParameters:
#     """Class storing Spectrum Dataset parameters

#     Parameters
#     ----------

#     """

#     def __init__(self,
#                  map_selection=None,
#                  methods=None,
#                  parameters=None,
#                  containment_correction=None,
#                  containment=None,
#                  use_region_center=None, 
#                  on_region_radius=None,
#                  acceptance=None,
#                  acceptance_off=None,
#                 ):
#         self.map_selection = map_selection
#         self.methods = methods
#         self.parameters = parameters
#         self.containment_correction = containment_correction
#         self.containment = containment
#         self.use_region_center = use_region_center
#         self.on_region_radius = on_region_radius
#         self.acceptance = acceptance
#         self.acceptance_off = acceptance_off


#     def __str__(self):
#         """Spectrum Dataset summary report (`str`)."""
#         ss = '*** Spectrum Dataset parameters summary ***\n'
#         ss += 'map_selection={}\n'.format(self.map_selection)
#         ss += 'methods={}\n'.format(self.methods)
#         ss += 'parameters={}\n'.format(self.parameters)
#         ss += 'containment_correction={}\n'.format(self.containment_correction)
#         ss += 'containment={}\n'.format(self.containment)
#         ss += 'use_region_center={}\n'.format(self.use_region_center)
#         ss += 'on_region_radius={}\n'.format(self.on_region_radius)
#         ss += 'acceptance={}\n'.format(self.acceptance)
#         ss += 'acceptance_off={}\n'.format(self.acceptance_off)
#         return ss
  

# def cut_energy_flux_points_datasets(datasets, e_ref_min=None, e_ref_max=None):
#     _datasets =  Datasets()
    
#     for dataset in datasets:
#         try: 
#             _dataset = cut_energy_flux_points_dataset(
#                 dataset=dataset, 
#                 e_ref_min=e_ref_min,
#                 e_ref_max=e_ref_max,
#             )
#             _datasets.append(_dataset)
#         except Exception as error:
#             print('\nUnable to cut {} FluxPointsDataset. An error has occurred: {}'.format(dataset.name, error))

 
#     return _datasets

# # @is_documented_by(Datasets)
# def write_datasets(datasets, filename=None, filename_models=None, overwrite=True):
#     """Write Datasets and Models to YAML file.

#         Parameters
#         ----------
#         overwrite : bool, optional
#             Overwrite existing file. Default is True.
#         """
    
#     if filename is None:
#         filename = "./datasets"
#     else: filename.mkdir(parents=True, exist_ok=True)
#     if filename_models:
#         filename_models.mkdir(parents=True, exist_ok=True)
        
#     datasets.write(filename=f"{filename}.yaml", filename_models=f"{filename_models}.yaml", overwrite=overwrite)

    
# # @is_documented_by([FluxPoints, FluxPointsDataset])
# def read_datasets(filename=None, filename_models=None):
#     """Read Datasets and Models from YAML file."""

#     if filename is None:
#         filename = "./datasets.yaml"
#     else: filename.mkdir(parents=True, exist_ok=True)
        
#     return Datasets.read(filename=filename, filename_models=filename_models)


# def create_spectrum_dataset_empty(geom, energy_axis_true, name="obs-0"):
#     """Create a MapDataset object with zero filled maps."""
#     return SpectrumDataset.create(
#         geom=geom, 
#         energy_axis_true=energy_axis_true,
#         name=name,
#     )


# def create_spectrum_dataset_onoff(dataset, acceptance, acceptance_off):
#     """Create a SpectrumDatasetOnOff from a `SpectrumDataset` dataset.""" 
#     return SpectrumDatasetOnOff.from_spectrum_dataset(
#         dataset=dataset, 
#         acceptance=acceptance, 
#         acceptance_off=acceptance_off,
#     )

# # In[ ]:
