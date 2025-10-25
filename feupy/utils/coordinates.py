# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities for coordinates."""

from astropy.coordinates import SkyCoord
import logging

log = logging.getLogger(__name__)

__all__ = [
    "convert_skycoord_to_dict",
    "convert_pos_config_to_skycoord",
#     "dict_to_skycoord",
]

def convert_skycoord_to_dict(position: SkyCoord) -> dict:
    """
    Convert a SkyCoord object to a dictionary with keys for longitude, latitude, and frame.

    Parameters
    ----------
    position : `~astropy.coordinates.SkyCoord`
        The position in the sky as a SkyCoord object.

    Returns
    -------
    dict
        Dictionary containing the following keys:
        - 'lon': Longitude of the position (`~astropy.coordinates.Angle`)
        - 'lat': Latitude of the position (`~astropy.coordinates.Angle`)
        - 'frame': Frame name of the coordinate system (str)
    
    Examples
    --------
    Convert a SkyCoord object to a dictionary:
    
    >>> from astropy.coordinates import SkyCoord
    >>> pos = SkyCoord(ra=10.684, dec=41.269, frame='icrs', unit='deg')
    >>> skycoord_to_dict(pos)
    {'lon': <Longitude 10.684 deg>, 'lat': <Latitude 41.269 deg>, 'frame': 'icrs'}
    """
    return {
        'lon': position.ra,
        'lat': position.dec,
        'frame': position.frame.name,
    }

def convert_pos_config_to_skycoord(pos_config) -> SkyCoord:
    """Convert a position configuration object to a SkyCoord object.

    Parameters:
    -----------
    pos_config : object
        An object containing position information with attributes:
        - lon : longitude or right ascension (RA) in astropy-compatible units
        - lat : latitude or declination (Dec) in astropy-compatible units
        - frame : coordinate frame (e.g., 'icrs', 'galactic')

    Returns:
    --------
    SkyCoord:
        A SkyCoord object with the specified longitude, latitude, and frame.

    Raises:
    -------
    AttributeError:
        If `pos_config` is missing required attributes.
    """
    try:
        return SkyCoord(pos_config.lon, pos_config.lat, frame=pos_config.frame)
    except AttributeError as e:
        log.error(f"Missing required attribute in pos_config: {e}")
        raise AttributeError(f"Invalid pos_config object: {e}")

# def skycoord_config_to_skycoord(pos_config):
#     return SkyCoord(pos_config.lon, pos_config.lat, frame=pos_config.frame)


# def dict_to_skycoord(pos_dict: dict):
#     return SkyCoord(pos_dict['lon'], pos_dict['lat'], frame=pos_dict['frame'])

# def skycoord_from_table(table):
#     keys = table.colnames

#     if {"RAJ2000", "DEJ2000"}.issubset(keys):
#         lon, lat, frame = "RAJ2000", "DEJ2000", "icrs"
#     elif {"RAJ2000", "DECJ2000"}.issubset(keys):
#         lon, lat, frame = "RAJ2000", "DECJ2000", "fk5"
#     elif {"RA", "DEC"}.issubset(keys):
#         lon, lat, frame = "RA", "DEC", "icrs"
#     elif {"ra", "dec"}.issubset(keys):
#         lon, lat, frame = "ra", "dec", "icrs"
#     else:
#         raise KeyError("No column GLON / GLAT or RA / DEC or RAJ2000 / DEJ2000 found.")

#     unit = table[lon].unit.to_string() if table[lon].unit else "deg"

#     return SkyCoord(table[lon], table[lat], unit=unit, frame=frame)