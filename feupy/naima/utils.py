# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Naima utilities for radiative models and spectral fitting.

This module provides helper functions and wrappers to use Naima models
within the Gammapy framework, including inverse Compton, pion decay models,
BIC calculation, and table conversions.

References
----------
- Naima: https://naima.readthedocs.io/
- Gammapy: https://docs.gammapy.org/
"""

import numpy as np
import astropy.units as u
from astropy.table import Table
from gammapy.estimators.map.core import DEFAULT_UNIT, OPTIONAL_QUANTITIES, REQUIRED_COLUMNS
from gammapy.modeling.models import NaimaSpectralModel, Models, SkyModel
from naima.plot import find_ML
import naima

REQUIRED_NAIMA_COLUMNS_NAMES = {
    'e_ref': 'energy',
    'e_min': 'energy_error_lo',
    'e_max': 'energy_error_hi',
    'dnde': 'flux',
    'dnde_err': 'flux_error',
    'dnde_errp': 'flux_error_hi',
    'dnde_errn': 'flux_error_lo',
    'dnde_ul': 'flux_ul',
    'e2dnde': 'flux',
    'e2dnde_err': 'flux_error',
    'e2dnde_errp': 'flux_error_hi',
    'e2dnde_errn': 'flux_error_lo',
    'e2dnde_ul': 'flux_ul',
    'is_ul': 'ul',
}

REQUIRED_NAIMA_COLUMNS = {
    'dnde': ['e_ref', 'dnde', 'dnde_err', 'dnde_errp', 'dnde_errn', 'dnde_ul', 'is_ul'],
    'e2dnde': ['e_ref', 'e2dnde', 'e2dnde_err', 'e2dnde_errp', 'e2dnde_errn', 'e2dnde_ul', 'is_ul'],
}

def get_inverse_compton_models(radiative_model, distance):
    """
    Create a list of SkyModels for inverse Compton emission components.

    Parameters
    ----------
    radiative_model : naima.radiative.RadiativeModel
        The Naima radiative model instance.
    distance : `~astropy.units.Quantity`
        Distance to the source.

    Returns
    -------
    models : `gammapy.modeling.models.Models`
        List of SkyModels for the total IC emission and individual seed photon fields.
    """
    models = Models()
    seeds = list(radiative_model.seed_photon_fields.keys())
    spectral_model = NaimaSpectralModel(radiative_model, distance=distance)
    model = SkyModel(spectral_model=spectral_model, name="IC (total)")
    models.append(model)

    for seed in seeds:
        spectral_model = NaimaSpectralModel(radiative_model, seed=seed, distance=distance)
        model = SkyModel(spectral_model=spectral_model, name=f"IC ({seed})")
        models.append(model)

    return models

def get_pion_decay_models(radiative_model, distance):
    """
    Create a SkyModel for pion decay emission.

    Parameters
    ----------
    radiative_model : naima.radiative.RadiativeModel
        The Naima radiative model instance.
    distance : `~astropy.units.Quantity`
        Distance to the source.

    Returns
    -------
    models : `gammapy.modeling.models.Models`
        List with a single SkyModel for pion decay emission.
    """
    models = Models()
    spectral_model = NaimaSpectralModel(radiative_model, distance=distance)
    model = SkyModel(spectral_model=spectral_model, name="Pion Decay")
    models.append(model)
    return models

def calc_BIC(sampler):
    """
    Calculate the Bayesian Information Criterion (BIC) from an MCMC sampler.

    Parameters
    ----------
    sampler : object
        Sampler object with `.data` and suitable for `naima.plot.find_ML`.

    Returns
    -------
    bic : float
        The Bayesian Information Criterion value.
    """
    MLp = find_ML(sampler, modelidx=0)[1]
    ML = find_ML(sampler, modelidx=0)[0]
    return len(MLp) * np.log(len(sampler.data)) - 2 * ML

def make_naima_tables(datasets, sed_type="dnde"):
    """
    Convert Gammapy datasets to Astropy tables formatted for Naima.

    Parameters
    ----------
    datasets : list
        List of Gammapy datasets containing spectral data.
    sed_type : str, optional
        SED type, either 'dnde' or 'e2dnde'. Default is 'dnde'.

    Returns
    -------
    tables : list of `~astropy.table.Table`
        List of tables formatted with Naima-compatible column names.
    """
    tables = []
    for dataset in datasets:
        table = Table()
        table.meta['name'] = dataset.name
        data = dataset.data.to_table(sed_type=sed_type)
        colnames = data.colnames
        columns = [x for x in REQUIRED_NAIMA_COLUMNS[sed_type] if x in colnames]
        columns_naima = [REQUIRED_NAIMA_COLUMNS_NAMES[x] for x in columns]
        for column, column_naima in zip(columns, columns_naima):
            table[column_naima] = data[column]
        tables.append(table)
    return tables

def compute_particle_distribution(radiative_model, E_min, E_max):
    """
    Compute the particle distribution and weighted particle distribution.

    Parameters
    ----------
    radiative_model : naima.radiative.RadiativeModel
        Radiative model providing the particle distribution.
    E_min : `~astropy.units.Quantity`
        Minimum particle energy.
    E_max : `~astropy.units.Quantity`
        Maximum particle energy.

    Returns
    -------
    energy : `~astropy.units.Quantity`
        Energy grid in TeV.
    particle_distribution : `~astropy.units.Quantity`
        Particle distribution evaluated on energy grid.
    particle_distribution_weighted : `~astropy.units.Quantity`
        Particle distribution weighted by E² and converted to erg.
    """
    energy = np.logspace(np.log10(E_min.to("TeV").value), np.log10(E_max.to("TeV").value), 100) * u.TeV
    particle_distribution = radiative_model.particle_distribution(energy)
    particle_distribution_weighted = (energy**2 * particle_distribution).to("erg")
    return energy, particle_distribution, particle_distribution_weighted

def get_sed_e2dnde(model_func, photon_energy, sed_type='e2dnde'):
    """
    Calculate the spectral energy distribution (SED) in E² dN/dE form.

    Parameters
    ----------
    model_func : callable
        Model function returning dN/dE.
    photon_energy : `~astropy.units.Quantity`
        Photon energy array.
    sed_type : str, optional
        Type of SED calculation (default 'e2dnde').

    Returns
    -------
    sed : `~astropy.units.Quantity`
        SED array in E² dN/dE units.
    """
    return (model_func * photon_energy) * photon_energy


# === Radiative Models ===

PionDecay_ECPL_p0 = np.array((46, 2.34, np.log10(80.0)))
PionDecay_ECPL_labels = ["log10(norm)", "index", "log10(cutoff)"]
proton_energy = np.logspace(-3, 2, 50) * u.TeV

def PionDecay_ECPL(pars, data):
    """
    Pion decay model with exponential cutoff power-law proton spectrum.

    Parameters
    ----------
    pars : array_like
        Model parameters: [log10(norm), index, log10(cutoff)].
    data : dict-like
        Data dictionary with energy information.

    Returns
    -------
    model : `~astropy.units.Quantity`
        Flux prediction.
    proton_distribution : tuple
        Tuple of (proton_energy, particle_distribution).
    Wp : `~astropy.units.Quantity`
        Total energy in protons above 1 TeV.
    """
    amplitude = 10 ** pars[0] / u.TeV
    alpha = pars[1]
    e_cutoff = 10 ** pars[2] * u.TeV
    ECPL = naima.models.ExponentialCutoffPowerLaw(amplitude, 30 * u.TeV, alpha, e_cutoff)
    PP = naima.models.PionDecay(ECPL, nh=1.0 * u.cm ** -3)
    model = PP.flux(data, distance=1.0 * u.kpc)
    proton_dist = PP.particle_distribution(proton_energy)
    Wp = PP.compute_Wp(Epmin=1 * u.TeV)
    return model, (proton_energy, proton_dist), Wp

def PionDecay_ECPL_lnprior(pars):
    """
    Log-prior for PionDecay_ECPL model.

    Parameters
    ----------
    pars : array_like
        Model parameters.

    Returns
    -------
    lnprior : float
        Logarithm of the prior probability.
    """
    return naima.uniform_prior(pars[1], -1, 5)


IC_We_p0 = np.array((40, 3.0, np.log10(30)))
IC_We_labels = ["log10(We)", "index", "log10(cutoff)"]

def IC_We(pars, data):
    """
    Inverse Compton model with electron energy cutoff and total energy We.

    Parameters
    ----------
    pars : array_like
        Model parameters: [log10(We), index, log10(cutoff)].
    data : dict-like
        Data dictionary with energy information.

    Returns
    -------
    model : `~astropy.units.Quantity`
        Flux prediction.
    electron_distribution : tuple
        Tuple of (electron_energy, particle_distribution).
    """
    We = 10 ** pars[0] * u.erg
    alpha = pars[1]
    e_cutoff = 10 ** pars[2] * u.TeV
    ECPL = naima.models.ExponentialCutoffPowerLaw(1 / u.eV, 10.0 * u.TeV, alpha, e_cutoff)
    IC = naima.models.InverseCompton(ECPL, seed_photon_fields=["CMB"])
    IC.set_We(We, Eemin=1 * u.TeV)
    model = IC.flux(data, distance=1.0 * u.kpc)
    elec_energy = np.logspace(11, 15, 100) * u.eV
    nelec = ECPL(elec_energy)
    return model, (elec_energy, nelec)

def IC_We_lnprior(pars):
    """
    Log-prior for IC_We model.

    Parameters
    ----------
    pars : array_like
        Model parameters.

    Returns
    -------
    lnprior : float
        Logarithm of the prior probability.
    """
    return naima.uniform_prior(pars[1], -1, 5)


# === Functional Models ===

ECPL_p0 = np.array((1e-12, 2.4, np.log10(15.0)))
ECPL_labels = ["norm", "index", "log10(cutoff)"]

def ECPL(pars, data):
    """
    Exponential cutoff power-law spectral model.

    Parameters
    ----------
    pars : array_like
        Model parameters: [norm, index, log10(cutoff)].
    data : dict-like
        Data dictionary with energy and flux.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Model flux prediction at given energies.
    """
    amplitude = pars[0] * data["flux"].unit
    alpha = pars[1]
    e_cutoff = (10 ** pars[2]) * u.TeV
    ECPL = naima.models.ExponentialCutoffPowerLaw(amplitude, 1 * u.TeV, alpha, e_cutoff)
    return ECPL(data)

def ECPL_lnprior(pars):
    """
    Log-prior for ECPL model.

    Parameters
    ----------
    pars : array_like
        Model parameters.

    Returns
    -------
    lnprior : float
        Logarithm of the prior probability.
    """
    return naima.uniform_prior(pars[0], 0.0, np.inf) + naima.uniform_prior(pars[1], -1, 5)


LP_p0 = np.array((1.5e-12, 2.7, 0.12))
LP_labels = ["norm", "alpha", "beta"]

def LP(pars, data):
    """
    Log-parabola spectral model.

    Parameters
    ----------
    pars : array_like
        Model parameters: [norm, alpha, beta].
    data : dict-like
        Data dictionary with energy and flux.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Model flux prediction at given energies.
    """
    amplitude = pars[0] * data["flux"].unit
    alpha = pars[1]
    beta = pars[2]
    LP = naima.models.LogParabola(amplitude, 1 * u.TeV, alpha, beta)
    return LP(data)

def LP_lnprior(pars):
    """
    Log-prior for LP model.

    Parameters
    ----------
    pars : array_like
        Model parameters.

    Returns
    -------
    lnprior : float
        Logarithm of the prior probability.
    """
    return naima.uniform_prior(pars[0], 0.0, np.inf) + naima.uniform_prior(pars[1], -1, 5)
