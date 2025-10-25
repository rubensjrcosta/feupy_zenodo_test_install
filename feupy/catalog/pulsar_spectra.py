# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities to fetch flux points tables and datasets for pulsar_spectra: A pulsar flux density catalogue and spectrum fitting repository (https://all-pulsar-spectra.readthedocs.io/en/latest/).

This module provides functions to:
- Read a pulsar spectra catalog.
- Generate flux points tables for individual pulsars.
- Group flux points tables by reference.
- Convert grouped tables into Gammapy datasets.
"""

import yaml
import logging
import numpy as np
from astropy.units import Quantity
from astropy.table import Table, Column
from gammapy.utils.scripts import make_path
from gammapy.datasets import Datasets, FluxPointsDataset
from gammapy.estimators import FluxPoints
from feupy.utils.units import Hz_to_eV, Jy_to_erg_by_cm2_s

# Set up logging
log = logging.getLogger(__name__)


__all__ = [
    "read_pulsar_catalog",
    "create_flux_points_table",
    "get_tables",
    "get_datasets",
]


def read_pulsar_catalog(filename="$FEUPY_DATA/data/catalogs/pulsar_spectra/pulsar_spectra.yaml"):
    """
    Read the pulsar spectra catalog from a YAML file.

    Parameters
    ----------
    filename : str
        Path to the pulsar spectra YAML file.

    Returns
    -------
    dict
        Dictionary containing the pulsar spectra catalog.
    """
    try:
        with open(make_path(filename), "r") as yaml_file:
            return yaml.safe_load(yaml_file)
    except yaml.YAMLError as error:
        log.error(f"YAML loading error: {error}")
        return {}


def create_flux_points_table(pulsar_jname):
    """
    Create a flux points table for a pulsar from the catalog.

    Parameters
    ----------
    pulsar_jname : str
        The pulsar J-name.

    Returns
    -------
    table : `~astropy.table.Table` or None
        Astropy table containing flux points data for the pulsar,
        or None if the pulsar is not found in the catalog.
    """
    catalog = read_pulsar_catalog()
    metadata = {
        "source_name": f"PSR {pulsar_jname}",
        "pulsar_jname": pulsar_jname,
        "catalog": "pulsar_spectra",
        "sed_type": "e2dnde",
        "comments": ["Reference: https://doi.org/10.1017/pasa.2022.52"],
        "bibcode": "2022PASA...39...56S",
    }

    try:
        # Extract pulsar data from the catalog
        freqs, bands, fluxs, flux_errs, refs = catalog[pulsar_jname]

        # Convert quantities to appropriate units
        freqs_mhz = Quantity(freqs, "MHz")
        fluxs_mjy = Quantity(fluxs, "mJy")
        flux_errs_mjy = Quantity(flux_errs, "mJy")

        # Create Astropy table with metadata
        table = Table(meta=metadata)
        table["ref"] = Column(data=np.array(refs, dtype="U20"), description="Reference label")
        table["e_ref"] = Column(data=Hz_to_eV(freqs_mhz), unit="eV", description="Reference energy", format=".3e")
        table["e2dnde"] = Column(
            data=Jy_to_erg_by_cm2_s(freqs_mhz, fluxs_mjy),
            unit="erg cm^-2 s^-1",
            description="Differential flux",
            format=".3e",
        )
        table["e2dnde_err"] = Column(
            data=Jy_to_erg_by_cm2_s(freqs_mhz, flux_errs_mjy),
            unit="erg cm^-2 s^-1",
            description="Differential flux uncertainty",
            format=".3e",
        )
        return table
    except KeyError:
        log.error(f"Error: Pulsar J-name '{pulsar_jname}' not found in the catalog.")
    except Exception as error:
        log.error(f"Error processing pulsar '{pulsar_jname}': {error}")
    return None


def get_flux_points_table(pulsar_jname):
    """
    Retrieve a flux points table for a specific pulsar J-name.

    Parameters
    ----------
    pulsar_jname : str
        The pulsar J-name.

    Returns
    -------
    table : `~astropy.table.Table` or None
        Flux points table for the pulsar, or None if the pulsar is not found.
    """
    return create_flux_points_table(pulsar_jname)


def get_tables(pulsar_jname):
    """
    Retrieve grouped flux points tables for a specific pulsar J-name.

    Parameters
    ----------
    pulsar_jname : str
        The pulsar J-name.

    Returns
    -------
    list of `~astropy.table.Table`
        List of grouped flux points tables for the pulsar,
        or an empty list if the pulsar is not found.
    """
    table = get_flux_points_table(pulsar_jname)
    if table is None:
        return []

    # Group table by reference
    grouped_tables = []
    grouped_data = table.group_by("ref")
    for group in grouped_data.groups:
        group_table = group[["e_ref", "e2dnde", "e2dnde_err"]]
        group_table.meta.update({"ref": group["ref"][0], **table.meta})
        grouped_tables.append(group_table)
    return grouped_tables


def get_datasets(pulsar_jname):
    """
    Retrieve datasets for a specific pulsar J-name.

    Parameters
    ----------
    pulsar_jname : str
        The pulsar J-name.

    Returns
    -------
    datasets : `~gammapy.datasets.Datasets`
        Collection of flux points datasets for the pulsar.
    """
    grouped_tables = get_tables(pulsar_jname)
    if not grouped_tables:
        return Datasets()

    # Create datasets from grouped tables
    datasets = Datasets()
    for group_table in grouped_tables:
        label = group_table.meta["ref"]
        flux_points = FluxPoints.from_table(table=group_table, sed_type=group_table.meta["sed_type"])
        datasets.append(FluxPointsDataset(data=flux_points, name=label))
    return datasets
