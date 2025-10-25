# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""2PC/3PC catalog and source classes."""

import astropy.units as u
from gammapy.estimators import FluxPoints
import logging

# Set up logging
log = logging.getLogger(__name__)

__all__ = [
    "get_flux_points_2PC",
    "get_flux_points_3PC",
]

def get_flux_points_2PC(source):
    """
    Generate flux points as `~gammapy.estimators.FluxPoints` for a 2PC source.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject2PC`
        Source object from the 2PC catalog.

    Returns
    -------
    flux_points : `~gammapy.estimators.FluxPoints`
        Flux points extracted from the source.
    """
    table = source.flux_points_table
    return FluxPoints.from_table(
        table=table,
        reference_model=source.spectral_model(),
    )

def get_flux_points_3PC(source, fit="auto"):
    """
    Generate flux points as `~gammapy.estimators.FluxPoints` for a 3PC source.

    In the 3PC, the Fermi-LAT collaboration attempted to fit a
    `~gammapy.modeling.models.SuperExpCutoffPowerLaw4FGLDR3SpectralModel` with the
    exponential index `index_2` either free or fixed to 2/3. These two models are referred
    to as "b free" and "b 23". For most pulsars, both models are available. However,
    in some cases, the "b free" model did not fit correctly.

    Parameters
    ----------
    source : `~gammapy.catalog.SourceCatalogObject3PC`
        Source object from the 3PC catalog.
    fit : str, optional
        Specifies which fitted model to return. Options are:
        - **"auto"** (default): Attempts to return the "b free" model first,
          falling back to "b 23" if "b free" is not available.
        - **"b free"**: Uses the "b free" model.
        - **"b 23"**: Uses the "b 23" model.

    Returns
    -------
    flux_points : `~gammapy.estimators.FluxPoints`
        Flux points extracted from the source.
    """
    table = source.flux_points_table
    return FluxPoints.from_table(
        table=table,
        reference_model=source.spectral_model(fit),
    )
