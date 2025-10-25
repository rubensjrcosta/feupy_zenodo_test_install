# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities for Geometry."""

# Imports
from astropy import units as u
from astropy.coordinates import SkyCoord
from gammapy.data import Observation, FixedPointingInfo, PointingMode
from gammapy.maps import MapAxis, RegionGeom
from regions import CircleSkyRegion


__all__ = [
    "create_energy_axis",
    "create_pointing",
    "create_pointing_position",
    "define_on_region",
    "create_region_geometry",
    # "GeometryParameters",
]


# Function Definitions

def create_energy_axis(energy_min, energy_max, nbin=5, per_decade=True, name="energy"):    
    """Create an energy axis for analysis."""
    return MapAxis.from_energy_bounds(
        energy_min=energy_min, 
        energy_max=energy_max, 
        nbin=nbin, 
        per_decade=per_decade, 
        name=name
    )


def create_pointing_position(position, position_angle, separation):
    """Calculate the pointing position based on a position, angle, and separation."""
    return position.directional_offset_by(position_angle, separation)


def create_pointing(pointing_position):
    """Create a pointing instance with a fixed position."""
    return FixedPointingInfo(
        mode=PointingMode.POINTING,
        fixed_icrs=pointing_position.icrs,
    )


def define_on_region(center, radius):
    """Define an on-region as a circular sky region."""
    return CircleSkyRegion(
        center=center, 
        radius=radius
    )


def create_region_geometry(on_region, axes):
    """Define the geometry for an analysis region."""
    return RegionGeom.create(
        region=on_region, 
        axes=axes
    )


# # Optional GeometryParameters Class
# class GeometryParameters:
#     """Container for geometry parameters.
# 
#     Parameters
#     ----------
#     e_reco_min : `~astropy.units.Quantity`
#         Minimal energy for simulation.
#     e_reco_max : `~astropy.units.Quantity`
#         Maximal energy for simulation.
#     nbin_reco : int
#         Number of bins for reconstructed energy.
#     e_true_min : `~astropy.units.Quantity`
#         Minimal true energy for simulation.
#     e_true_max : `~astropy.units.Quantity`
#         Maximal true energy for simulation.
#     nbin_true : int
#         Number of bins for true energy.
#     """
#     @u.quantity_input(
#         e_reco_min=u.eV, 
#         e_reco_max=u.eV,
#         e_true_min=u.eV, 
#         e_true_max=u.eV
#     )
#     def __init__(self,
#                  e_reco_min=None,
#                  e_reco_max=None,
#                  nbin_reco: int=None,
#                  e_true_min=None,
#                  e_true_max=None,
#                  nbin_true: int=None,
#                 ):
#         self.e_reco_min = Quantity(e_reco_min, "TeV")
#         self.e_reco_max = Quantity(e_reco_max, "TeV")
#         self.nbin_reco = nbin_reco
#         self.e_true_min = Quantity(e_true_min, "TeV")
#         self.e_true_max = Quantity(e_true_max, "TeV")
#         self.nbin_true = nbin_true
# 
#     def __str__(self):
#         """Return a summary of geometry parameters."""
#         ss = '*** Basic parameters ***\n\n'
#         ss += 'e_reco_min = {:.2f}\n'.format(self.e_reco_min).replace(' ', '')
#         ss += 'e_reco_max = {:.2f}\n'.format(self.e_reco_max).replace(' ', '')
#         ss += 'nbin_reco = {}\n'.format(self.nbin_reco)
#         ss += 'e_true_min = {:.2f}\n'.format(self.e_true_min).replace(' ', '')
#         ss += 'e_true_max = {:.2f}\n'.format(self.e_true_max).replace(' ', '')
#         ss += 'nbin_true = {}\n'.format(self.nbin_true)
#         return ss
