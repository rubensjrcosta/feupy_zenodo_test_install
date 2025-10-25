# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" Utility functions and constants for unit conversions in Feupy.

You have to import sub-modules of `feupy.utils` directly,
the `feupy.utils` namespace is empty.

This module provides a collection of utility functions and constants that are used throughout Feupy for handling unit conversions, particularly in relation to spectral energy distributions (SEDs) and flux calculations. The functions included help with converting between commonly used units such as frequency (Hz), energy (eV), and flux density (Jy) to physical flux units (erg/cm^2/s).

Key functions include:

Jy_to_erg_by_cm2_s: Converts flux density from Jansky (Jy) to physical flux in erg/cm^2/s.

Hz_to_eV: Converts frequency (Hz) to energy (eV).

xaxis_units_to_energy_label: Generates an appropriate energy label for plotting based on the provided X-axis units.

Additionally, the module defines constants for default axis labels used in plotting Spectral Energy Distributions (SEDs).

Examples:

from gammapy.utils.units import Jy_to_erg_by_cm2_s, Hz_to_eV

freq = 1e9 * u.Hz  # Frequency in Hz
flux = 1 * u.mJy  # Flux in mJy
energy_flux = Jy_to_erg_by_cm2_s(freq, flux)

freq_to_energy = Hz_to_eV(freq)
"""