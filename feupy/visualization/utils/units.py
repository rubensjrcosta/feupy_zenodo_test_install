# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Units-related utilities and constants for Feupy."""

from astropy.units import Quantity
from astropy import units as u
from gammapy.maps.axes import UNIT_STRING_FORMAT
from gammapy.estimators.map.core import DEFAULT_UNIT

# Define default Y-axis labels for plotting SEDs
DEFAULT_YAXIS_LABEL = {
    'e2dnde': f"[{DEFAULT_UNIT['e2dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {E^{2}\\,\Phi(E)\\, ['),
    'dnde': f"[{DEFAULT_UNIT['dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {\Phi(E)\\, [')
}

# # Define default Y-axis labels for energy flux SEDs
# YAXIS_LABEL = {
#     'e2dnde': f"[{DEFAULT_UNIT['e2dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {E^{2}\\,J(E)\\, ['),
#     'dnde': f"[{DEFAULT_UNIT['dnde'].to_string(UNIT_STRING_FORMAT)}]".replace('[$\\mathrm{', '$\\rm {J(E)\\, [')
# }

# Define default X-axis labels for energy axes
DEFAULT_XAXIS_LABEL = {
    'erg': f"Energy [{u.Unit('erg').to_string(UNIT_STRING_FORMAT)}]",
    'TeV': f"Energy [{u.Unit('TeV').to_string(UNIT_STRING_FORMAT)}]"
}
