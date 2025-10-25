# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Statistics."""
from scipy.stats import chi2, norm
from gammapy.stats import WStatCountsStatistic
from feupy.utils.scripts import is_documented_by
import numpy as np
from gammapy.modeling import Fit
from gammapy.datasets import Datasets, FluxPointsDataset
from gammapy.estimators import FluxPoints
from gammapy.modeling.models import SkyModel

__all__ = [
    "StatisticalUtilityFunctions",
    "calculate_AIC",
    "calculate_relative_AIC",
    "fit_spectral_model_to_flux_points",
]


class StatisticalUtilityFunctions:
    """Statistical utility functions
    
    See: https://docs.gammapy.org/1.1/user-guide/stats/
    
    StatisticalUtilityFunctions is represented by `~feupy.analysis.simulation.stats`.
    """    
    
    def __init__(self):
        pass
        
#         self.irfs_list = None
#         self.irfs_label_list = None
    
    @classmethod
    def sigma_to_ts(cls, sigma, df=1):
        """Convert sigma to delta ts"""
        p_value = 2 * norm.sf(sigma)
        return chi2.isf(p_value, df=df)
    
    @classmethod
    def ts_to_sigma(cls, ts, df=1):
        """Convert delta ts to sigma"""
        p_value = chi2.sf(ts, df=df)
        return norm.isf(0.5 * p_value)
    
    @classmethod
    def calculate_sensitivity_lima(cls, n_on_events, n_background, alpha, n_bins_energy,
                               n_bins_gammaness, n_bins_theta2):
        """
        Sensitivity calculation using the Li & Ma formula
        eq. 17 of Li & Ma (1983).
        https://ui.adsabs.harvard.edu/abs/1983ApJ...272..317L/abstract

        We calculate the sensitivity in bins of energy, gammaness and
        theta2

        Parameters
        ---------
        n_on_events:   `numpy.ndarray` number of ON events in the signal region
        n_background:   `numpy.ndarray` number of events in the background region
        alpha: `float` inverse of the number of off positions
        n_bins_energy: `int` number of bins in energy
        n_bins_gammaness: `int` number of bins in gammaness
        n_bins_theta2: `int` number of bins in theta2

        Returns
        ---------
        sensitivity: `numpy.ndarray` sensitivity in percentage of Crab units
        n_excesses_5sigma: `numpy.ndarray` number of excesses corresponding to 
                    a 5 sigma significance

        """

        stat = WStatCountsStatistic(n_on=n_on_events,
                                    n_off=n_background,
                                    alpha=alpha)

        n_excesses_5sigma = stat.excess_matching_significance(5)

        for i in range(0, n_bins_energy):
            for j in range(0, n_bins_gammaness):
                for k in range(0, n_bins_theta2):
                    if n_excesses_5sigma[i][j][k] < 10:
                        n_excesses_5sigma[i][j][k] = 10

                    if n_excesses_5sigma[i, j,
                                         k] < 0.05 * n_background[i][j][k] / 5:
                        n_excesses_5sigma[i, j,
                                          k] = 0.05 * n_background[i][j][k] / 5

        sensitivity = n_excesses_5sigma / n_on_events * 100  # percentage of Crab

        return n_excesses_5sigma, sensitivity

    
    @classmethod
    def calculate_sensitivity_lima_ebin(cls, n_excesses, n_background, alpha, n_bins_energy):
        """
        Sensitivity calculation using the Li & Ma formula
        eq. 17 of Li & Ma (1983).
        https://ui.adsabs.harvard.edu/abs/1983ApJ...272..317L/abstract

        Parameters
        ---------
        n_excesses:   `numpy.ndarray` number of excess events in the signal region
        n_background: `numpy.ndarray` number of events in the background region
        alpha:        `float` inverse of the number of off positions
        n_bins_energy:`int` number of bins in energy

        Returns
        ---------
        sensitivity: `numpy.ndarray` sensitivity in percentage of Crab units
        n_excesses_5sigma: `numpy.ndarray` number of excesses corresponding to 
                    a 5 sigma significance

        """

        if any(len(a) != n_bins_energy for a in (n_excesses, n_background, alpha)):
            raise ValueError(
                'Excess, background and alpha arrays must have the same length')

        stat = WStatCountsStatistic(
            n_on=np.ones_like(n_background),
            n_off=n_background,
            alpha=alpha)

        n_excesses_5sigma = stat.n_sig_matching_significance(5)

        for i in range(0, n_bins_energy):
            # If the excess needed to get 5 sigma is less than 10,
            # we force it to be at least 10
            if n_excesses_5sigma[i] < 10:
                n_excesses_5sigma[i] = 10
            # If the excess needed to get 5 sigma is less than 5%
            # of the background, we force it to be at least 5% of
            # the background
            if n_excesses_5sigma[i] < 0.05 * n_background[i] * alpha[i]:
                n_excesses_5sigma[i] = 0.05 * n_background[i] * alpha[i]

        sensitivity = n_excesses_5sigma / n_excesses * 100  # percentage of Crab

        return n_excesses_5sigma, sensitivity
    
    @classmethod
    def compute_significance(cls, dataset_onoff, alpha=0.2):
        # Class to compute statistics for Poisson distributed variable with unknown background.
        return WStatCountsStatistic(
            n_on=sum(dataset_onoff.counts.data), 
            n_off=sum(dataset_onoff.counts_off.data), 
            alpha=alpha).sqrt_ts
    
    @classmethod
    @is_documented_by(WStatCountsStatistic)
    def compute_wstat(cls, dataset_onoff, alpha=0.2):
        # Class to compute statistics for Poisson distributed variable with unknown background.
        return WStatCountsStatistic(
            n_on=sum(dataset_onoff.counts.data), 
            n_off=sum(dataset_onoff.counts_off.data), 
            alpha=alpha)


def calculate_relative_AIC(datasets, fit_result_H0, fit_result_H1):
    """
    Calculate the relative difference in Akaike Information Criterion (AIC).

    Parameters
    ----------
    datasets : `Datasets` or list of `Dataset`
        Datasets used in the fit.
    fit_result_H0 : `FitResult`
        Fit result of the hypothesis 0.
    fit_result_H1 : `FitResult`
        Fit result of the hypothesis 1.
        
    Returns
    -------
    delta_AIC : `float`
        Relative difference in AIC results
    """

    AICc_0 = calculate_AIC(datasets, fit_result_H0)
    AICc_1 = calculate_AIC(datasets, fit_result_H1)


    # Calculate the relative difference in AIC
    delta_AIC = (1 - AICc_1 / AICc_0) * 100
    print(f"Delta AIC_{tag_H1} = {delta_AIC:.2f}%")

    return delta_AIC

def calculate_AIC(datasets, fit_result):
    """
    Calculate Akaike Information Criterion (AIC).

    Parameters
    ----------
    datasets : `Datasets` or list of `Dataset`
        Datasets used in the fit.
    fit_result : `FitResult`
        Fit result.        
    Returns
    -------
    AIC: `float`
        Akaike information criterion (AIC) result.
    """
    # Check if both fits were successful
    if not fit_result.success:
        raise ValueError("Fit for Hypothesis was not successful.")

    # Extract model tags for labeling
    tag_H0 = fit_result.models[0].spectral_model.tag[1].upper()

    # Calculate the number of data points (no upper limits)
    N_pt = sum(
        int(np.sum(~dataset.data.is_ul.data))  # Ensure scalar integer
        for dataset in datasets
    )

    # Extract Wstat values
    Wstat_0 = float(fit_result.total_stat)  # Ensure scalar

    # Extract number of free parameters
    k_0 = int(len(fit_result.models.parameters.free_parameters.names))

    # Calculate AIC and corrected AIC
    AIC_0 = Wstat_0 + 2 * k_0

    AICc_0 = AIC_0 + ((2 * k_0**2 + 2 * k_0) / (N_pt - k_0 - 1))

    print(f"AIC_{tag_H0} = {AICc_0:.2f}")


    return AICc_0

def fit_spectral_model_to_flux_points(flux_points_table, spectral_model):
    """Fit a spectral model to flux points.

    This function takes a table of flux points and a spectral model, fits the model to
    the flux points using a likelihood approach, and returns the fitted spectral model.

    Parameters
    ----------
    flux_points_table : `~astropy.table.Table`
        Table containing the flux points data to be fitted. The table must have
        columns that are compatible with `gammapy.estimators.FluxPoints`.
    
    spectral_model : `~gammapy.modeling.models.SpectralModel`
        The spectral model to be fitted to the flux points. This model should be
        initialized with starting parameters.

    Returns
    -------
    spectral_model : `~gammapy.modeling.models.SpectralModel`
        The spectral model with optimized parameters after fitting to the flux points.

    Examples
    --------
    Fit a log-parabola spectral model to flux points from a table:

    >>> from astropy.table import Table
    >>> from gammapy.modeling.models import LogParabolaSpectralModel
    >>> flux_points_table = Table.read("flux_points.ecsv")
    >>> spectral_model = LogParabolaSpectralModel()
    >>> fitted_model = fit_spectral_model_to_flux_points(flux_points_table, spectral_model)
    >>> print(fitted_model)
    """
    flux_points = FluxPoints.from_table(flux_points_table)
    datasets = Datasets(FluxPointsDataset(data=flux_points))
    model = SkyModel(spectral_model=spectral_model)
    datasets.models = model
    fitter = Fit()
    result = fitter.run(datasets=datasets)

    return model.spectral_model

