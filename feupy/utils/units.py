# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Units-related utilities and constants for Feupy."""

from astropy.units import Quantity
from astropy import units as u

UNIT_DEG = 'deg' 
# Degree unit
FRAME_ICRS = 'icrs'
# Frame icrs
FRAME_FK5 = 'fk5'
# Frame fk5

    
ENERGY_COLUMNS = {
    "dnde": (["e_ref"], ),
    "e2dnde": ["e_ref"],
    "flux": ["e_min", "e_max", "flux_sensitivity", "flux_err", "flux_errp", "flux_errn"],
    "eflux": ["e_min", "e_max", "eflux_sensitivity", "eflux_err", "eflux_errp", "eflux_errn"],
    # Extended required columns
    "likelihood": [
        "e_min",
        "e_max",
        "likelihood_ratio",
        "delta_log_likelihood"
    ],
}

SED_COLUMNS = {
    "dnde": ["dnde", "dnde_err", "dnde_errp", "dnde_errn", "dnde_ul"],
    "e2dnde": ["e2dnde", "e2dnde_err", "e2dnde_errp", "e2dnde_errn", "e2dnde_ul"],
    "flux": ["flux", "flux_err", "flux_errp", "flux_errn", "flux_ul", "flux_sensitivity"],
    "eflux": ["eflux", "eflux_err", "eflux_errp", "eflux_errn", "eflux_ul", "eflux_sensitivity"],
    # Extended required columns
    "likelihood": [
        "e_min",
        "e_max",
        "e_ref",
        "ref_dnde",
        "ref_flux",
        "ref_eflux",
        "norm",
        "norm_err", "norm_errn", "norm_errp", "norm_ul",
        "likelihood_ratio",
        "delta_log_likelihood"
    ],
}

DEFAULT_ENERGY_UNIT = {
    "dnde": u.Unit("TeV"),
    "e2dnde": u.Unit("TeV"),
    "flux": u.Unit("TeV"),
    "eflux": u.Unit("TeV"),
    "norm": u.Unit(""),
    # Possible extensions
    "eflux_density": u.Unit("erg cm-2 s-1 TeV-1"),
    "counts": u.Unit(""),
}

DEFAULT_SED_UNIT = {
    "dnde": u.Unit("cm-2 s-1 TeV-1"),
    "e2dnde": u.Unit("erg cm-2 s-1"),
    "flux": u.Unit("cm-2 s-1"),
    "eflux": u.Unit("erg cm-2 s-1"),
    "norm": u.Unit(""),
    # Possible extensions
    "eflux_density": u.Unit("erg cm-2 s-1 TeV-1"),
    "counts": u.Unit("cm-2 s-1"),
}

def Jy_to_erg_by_cm2_s(freq, flux):
    """
    Convert flux density from Jansky (Jy) to flux in erg/cm^2/s.

    Parameters
    ----------
    freq : `~astropy.units.Quantity`
        Frequency in units of Hertz (Hz).
    flux : `~astropy.units.Quantity`
        Flux density in units of milliJansky (mJy).
        
    Returns
    -------
    flux : `~astropy.units.Quantity`
        Flux in erg/cm^2/s.
    """
    return Quantity(flux.to(u.mJy) * freq.to(u.Hz), 'Jy Hz').to('erg cm-2 s-1')


def Hz_to_eV(freq):
    """
    Convert frequency (Hz) to energy (eV).

    Parameters
    ----------
    freq : `~astropy.units.Quantity`
        Frequency in units of Hertz (Hz).
        
    Returns
    -------
    energy : `~astropy.units.Quantity`
        Energy in units of electron volts (eV).
    """
    return Quantity(freq).to(u.eV, equivalencies=u.spectral())


# from astropy.units import Quantity
# from gammapy.maps.axes import UNIT_STRING_FORMAT
# from gammapy.estimators.map.core import DEFAULT_UNIT





# # Define default Y-axis labels for plotting SEDs
# DEFAULT_YAXIS_LABEL = {
#     'e2dnde': f"[{DEFAULT_UNIT['e2dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {E^{2}\\,\Phi(E)\\, ['),
#     'dnde': f"[{DEFAULT_UNIT['dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {\Phi(E)\\, [')
# }

# # # Define default Y-axis labels for energy flux SEDs
# # YAXIS_LABEL = {
# #     'e2dnde': f"[{DEFAULT_UNIT['e2dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {E^{2}\\,J(E)\\, ['),
# #     'dnde': f"[{DEFAULT_UNIT['dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {J(E)\\, [')
# # }

# # Define default X-axis labels for energy axes
# DEFAULT_XAXIS_LABEL = {
#     'erg': f"Energy [{u.Unit('erg').to_string(UNIT_STRING_FORMAT)}]",
#     'TeV': f"Energy [{u.Unit('TeV').to_string(UNIT_STRING_FORMAT)}]"
# }


# def xaxis_units_to_energy_label(xaxis_units):
#     """
#     Generate an energy label from X-axis units for plotting.

#     Parameters
#     ----------
#     xaxis_units : `~astropy.units.Unit`
#         The unit for the X-axis, typically an energy unit.
        
#     Returns
#     -------
#     str
#         The label for the X-axis in the appropriate unit format.
        
#     Notes
#     -----
#     If the provided unit is not compatible with energy units, an error will be printed.
#     """
#     try:
#         (1 * xaxis_units).to(u.eV)
#         return f"Energy [{xaxis_units.to_string(UNIT_STRING_FORMAT)}]"
#     except Exception as error:
#         print(f"Error: {error}")
