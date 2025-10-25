# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""CTAO IRFs class."""

from gammapy.irf import load_irf_dict_from_file
from gammapy.data import observatory_locations
from astropy.units import Quantity
from astropy.coordinates import Angle
import astropy.units as u
from gammapy.data import observatory_locations
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import numpy as np
from datetime import datetime, timedelta


__all__ = ["Irfs"]

# Define options for IRFs
_IRFS_OPTIONS = [
    ['South', 'AverageAz', '20deg', '0.5h'],
     ['South', 'AverageAz', '20deg', '5h'],
     ['South', 'AverageAz', '20deg', '50h'],
     ['South', 'NorthAz', '20deg', '0.5h'],
     ['South', 'NorthAz', '20deg', '5h'],
     ['South', 'NorthAz', '20deg', '50h'],
     ['South', 'SouthAz', '20deg', '0.5h'],
     ['South', 'SouthAz', '20deg', '5h'],
     ['South', 'SouthAz', '20deg', '50h'],
     ['South', 'AverageAz', '40deg', '0.5h'],
     ['South', 'AverageAz', '40deg', '5h'],
     ['South', 'AverageAz', '40deg', '50h'],
     ['South', 'NorthAz', '40deg', '0.5h'],
     ['South', 'NorthAz', '40deg', '5h'],
     ['South', 'NorthAz', '40deg', '50h'],
     ['South', 'SouthAz', '40deg', '0.5h'],
     ['South', 'SouthAz', '40deg', '5h'],
     ['South', 'SouthAz', '40deg', '50h'],
     ['South', 'AverageAz', '60deg', '0.5h'],
     ['South', 'AverageAz', '60deg', '5h'],
     ['South', 'AverageAz', '60deg', '50h'],
     ['South', 'NorthAz', '60deg', '0.5h'],
     ['South', 'NorthAz', '60deg', '5h'],
     ['South', 'NorthAz', '60deg', '50h'],
     ['South', 'SouthAz', '60deg', '0.5h'],
     ['South', 'SouthAz', '60deg', '5h'],
     ['South', 'SouthAz', '60deg', '50h'],
     ['South-SSTSubArray', 'AverageAz', '20deg', '0.5h'],
     ['South-SSTSubArray', 'AverageAz', '20deg', '5h'],
     ['South-SSTSubArray', 'AverageAz', '20deg', '50h'],
     ['South-SSTSubArray', 'NorthAz', '20deg', '0.5h'],
     ['South-SSTSubArray', 'NorthAz', '20deg', '5h'],
     ['South-SSTSubArray', 'NorthAz', '20deg', '50h'],
     ['South-SSTSubArray', 'SouthAz', '20deg', '0.5h'],
     ['South-SSTSubArray', 'SouthAz', '20deg', '5h'],
     ['South-SSTSubArray', 'SouthAz', '20deg', '50h'],
     ['South-SSTSubArray', 'AverageAz', '40deg', '0.5h'],
     ['South-SSTSubArray', 'AverageAz', '40deg', '5h'],
     ['South-SSTSubArray', 'AverageAz', '40deg', '50h'],
     ['South-SSTSubArray', 'NorthAz', '40deg', '0.5h'],
     ['South-SSTSubArray', 'NorthAz', '40deg', '5h'],
     ['South-SSTSubArray', 'NorthAz', '40deg', '50h'],
     ['South-SSTSubArray', 'SouthAz', '40deg', '0.5h'],
     ['South-SSTSubArray', 'SouthAz', '40deg', '5h'],
     ['South-SSTSubArray', 'SouthAz', '40deg', '50h'],
     ['South-SSTSubArray', 'AverageAz', '60deg', '0.5h'],
     ['South-SSTSubArray', 'AverageAz', '60deg', '5h'],
     ['South-SSTSubArray', 'AverageAz', '60deg', '50h'],
     ['South-SSTSubArray', 'NorthAz', '60deg', '0.5h'],
     ['South-SSTSubArray', 'NorthAz', '60deg', '5h'],
     ['South-SSTSubArray', 'NorthAz', '60deg', '50h'],
     ['South-SSTSubArray', 'SouthAz', '60deg', '0.5h'],
     ['South-SSTSubArray', 'SouthAz', '60deg', '5h'],
     ['South-SSTSubArray', 'SouthAz', '60deg', '50h'],
     ['South-MSTSubArray', 'AverageAz', '20deg', '0.5h'],
     ['South-MSTSubArray', 'AverageAz', '20deg', '5h'],
     ['South-MSTSubArray', 'AverageAz', '20deg', '50h'],
     ['South-MSTSubArray', 'NorthAz', '20deg', '0.5h'],
     ['South-MSTSubArray', 'NorthAz', '20deg', '5h'],
     ['South-MSTSubArray', 'NorthAz', '20deg', '50h'],
     ['South-MSTSubArray', 'SouthAz', '20deg', '0.5h'],
     ['South-MSTSubArray', 'SouthAz', '20deg', '5h'],
     ['South-MSTSubArray', 'SouthAz', '20deg', '50h'],
     ['South-MSTSubArray', 'AverageAz', '40deg', '0.5h'],
     ['South-MSTSubArray', 'AverageAz', '40deg', '5h'],
     ['South-MSTSubArray', 'AverageAz', '40deg', '50h'],
     ['South-MSTSubArray', 'NorthAz', '40deg', '0.5h'],
     ['South-MSTSubArray', 'NorthAz', '40deg', '5h'],
     ['South-MSTSubArray', 'NorthAz', '40deg', '50h'],
     ['South-MSTSubArray', 'SouthAz', '40deg', '0.5h'],
     ['South-MSTSubArray', 'SouthAz', '40deg', '5h'],
     ['South-MSTSubArray', 'SouthAz', '40deg', '50h'],
     ['South-MSTSubArray', 'AverageAz', '60deg', '0.5h'],
     ['South-MSTSubArray', 'AverageAz', '60deg', '5h'],
     ['South-MSTSubArray', 'AverageAz', '60deg', '50h'],
     ['South-MSTSubArray', 'NorthAz', '60deg', '0.5h'],
     ['South-MSTSubArray', 'NorthAz', '60deg', '5h'],
     ['South-MSTSubArray', 'NorthAz', '60deg', '50h'],
     ['South-MSTSubArray', 'SouthAz', '60deg', '0.5h'],
     ['South-MSTSubArray', 'SouthAz', '60deg', '5h'],
     ['South-MSTSubArray', 'SouthAz', '60deg', '50h'],
     ['North', 'AverageAz', '20deg', '0.5h'],
     ['North', 'AverageAz', '20deg', '5h'],
     ['North', 'AverageAz', '20deg', '50h'],
     ['North', 'NorthAz', '20deg', '0.5h'],
     ['North', 'NorthAz', '20deg', '5h'],
     ['North', 'NorthAz', '20deg', '50h'],
     ['North', 'SouthAz', '20deg', '0.5h'],
     ['North', 'SouthAz', '20deg', '5h'],
     ['North', 'SouthAz', '20deg', '50h'],
     ['North', 'AverageAz', '40deg', '0.5h'],
     ['North', 'AverageAz', '40deg', '5h'],
     ['North', 'AverageAz', '40deg', '50h'],
     ['North', 'NorthAz', '40deg', '0.5h'],
     ['North', 'NorthAz', '40deg', '5h'],
     ['North', 'NorthAz', '40deg', '50h'],
     ['North', 'SouthAz', '40deg', '0.5h'],
     ['North', 'SouthAz', '40deg', '5h'],
     ['North', 'SouthAz', '40deg', '50h'],
     ['North', 'AverageAz', '60deg', '0.5h'],
     ['North', 'AverageAz', '60deg', '5h'],
     ['North', 'AverageAz', '60deg', '50h'],
     ['North', 'NorthAz', '60deg', '0.5h'],
     ['North', 'NorthAz', '60deg', '5h'],
     ['North', 'NorthAz', '60deg', '50h'],
     ['North', 'SouthAz', '60deg', '0.5h'],
     ['North', 'SouthAz', '60deg', '5h'],
     ['North', 'SouthAz', '60deg', '50h'],
     ['North-MSTSubArray', 'AverageAz', '20deg', '0.5h'],
     ['North-MSTSubArray', 'AverageAz', '20deg', '5h'],
     ['North-MSTSubArray', 'AverageAz', '20deg', '50h'],
     ['North-MSTSubArray', 'NorthAz', '20deg', '0.5h'],
     ['North-MSTSubArray', 'NorthAz', '20deg', '5h'],
     ['North-MSTSubArray', 'NorthAz', '20deg', '50h'],
     ['North-MSTSubArray', 'SouthAz', '20deg', '0.5h'],
     ['North-MSTSubArray', 'SouthAz', '20deg', '5h'],
     ['North-MSTSubArray', 'SouthAz', '20deg', '50h'],
     ['North-MSTSubArray', 'AverageAz', '40deg', '0.5h'],
     ['North-MSTSubArray', 'AverageAz', '40deg', '5h'],
     ['North-MSTSubArray', 'AverageAz', '40deg', '50h'],
     ['North-MSTSubArray', 'NorthAz', '40deg', '0.5h'],
     ['North-MSTSubArray', 'NorthAz', '40deg', '5h'],
     ['North-MSTSubArray', 'NorthAz', '40deg', '50h'],
     ['North-MSTSubArray', 'SouthAz', '40deg', '0.5h'],
     ['North-MSTSubArray', 'SouthAz', '40deg', '5h'],
     ['North-MSTSubArray', 'SouthAz', '40deg', '50h'],
     ['North-MSTSubArray', 'AverageAz', '60deg', '0.5h'],
     ['North-MSTSubArray', 'AverageAz', '60deg', '5h'],
     ['North-MSTSubArray', 'AverageAz', '60deg', '50h'],
     ['North-MSTSubArray', 'NorthAz', '60deg', '0.5h'],
     ['North-MSTSubArray', 'NorthAz', '60deg', '5h'],
     ['North-MSTSubArray', 'NorthAz', '60deg', '50h'],
     ['North-MSTSubArray', 'SouthAz', '60deg', '0.5h'],
     ['North-MSTSubArray', 'SouthAz', '60deg', '5h'],
     ['North-MSTSubArray', 'SouthAz', '60deg', '50h'],
     ['North-LSTSubArray', 'AverageAz', '20deg', '0.5h'],
     ['North-LSTSubArray', 'AverageAz', '20deg', '5h'],
     ['North-LSTSubArray', 'AverageAz', '20deg', '50h'],
     ['North-LSTSubArray', 'NorthAz', '20deg', '0.5h'],
     ['North-LSTSubArray', 'NorthAz', '20deg', '5h'],
     ['North-LSTSubArray', 'NorthAz', '20deg', '50h'],
     ['North-LSTSubArray', 'SouthAz', '20deg', '0.5h'],
     ['North-LSTSubArray', 'SouthAz', '20deg', '5h'],
     ['North-LSTSubArray', 'SouthAz', '20deg', '50h'],
     ['North-LSTSubArray', 'AverageAz', '40deg', '0.5h'],
     ['North-LSTSubArray', 'AverageAz', '40deg', '5h'],
     ['North-LSTSubArray', 'AverageAz', '40deg', '50h'],
     ['North-LSTSubArray', 'NorthAz', '40deg', '0.5h'],
     ['North-LSTSubArray', 'NorthAz', '40deg', '5h'],
     ['North-LSTSubArray', 'NorthAz', '40deg', '50h'],
     ['North-LSTSubArray', 'SouthAz', '40deg', '0.5h'],
     ['North-LSTSubArray', 'SouthAz', '40deg', '5h'],
     ['North-LSTSubArray', 'SouthAz', '40deg', '50h'],
     ['North-LSTSubArray', 'AverageAz', '60deg', '0.5h'],
     ['North-LSTSubArray', 'AverageAz', '60deg', '5h'],
     ['North-LSTSubArray', 'AverageAz', '60deg', '50h'],
     ['North-LSTSubArray', 'NorthAz', '60deg', '0.5h'],
     ['North-LSTSubArray', 'NorthAz', '60deg', '5h'],
     ['North-LSTSubArray', 'NorthAz', '60deg', '50h'],
     ['North-LSTSubArray', 'SouthAz', '60deg', '0.5h'],
     ['North-LSTSubArray', 'SouthAz', '60deg', '5h'],
     ['North-LSTSubArray', 'SouthAz', '60deg', '50h']
]

class Irfs:
    """Class to handle Instrument Response Functions (IRFs) for CTAO."""
    
    IRFS_OPTIONS = _IRFS_OPTIONS
    IRF_VERSION = "prod5 v0.1"
    
    _SITE_ARRAY = {
        'South': '14MSTs37SSTs', 
        'South-SSTSubArray': '37SSTs', 
        'South-MSTSubArray': '14MSTs', 
        'North': '4LSTs09MSTs', 
        'North-MSTSubArray': '09MSTs',
        'North-LSTSubArray': '4LSTs'
    }
    _OBS_TIME = {'0.5h': '1800s', '5h': '18000s', '50h': '180000s'}
    
    _DIR_FITS = '$FEUPY_DATA/data/irfs/cta-prod5-zenodo-v0.1/fits/'

    def __init__(self):
        self.irfs = None
        self.irfs_label = None
        self.obs_loc = None

    @classmethod
    def get_irfs(cls, irfs_opt):
        """Load IRF file based on specified options."""
        dir_fits = f'CTA-Performance-prod5-v0.1-{irfs_opt[0]}-{irfs_opt[2]}.FITS/'
        isite = irfs_opt[0].split('-')[0]
        irfs_file_name = (
            f'Prod5-{isite}-{irfs_opt[2]}-{irfs_opt[1]}-'
            f'{cls._SITE_ARRAY[irfs_opt[0]]}.{cls._OBS_TIME[irfs_opt[3]]}-v0.1.fits.gz'
        )
        file_path = f'{cls._DIR_FITS}{dir_fits}{irfs_file_name}'
        cls.irfs = load_irf_dict_from_file(file_path)
        cls.irfs_label = cls.get_irfs_label(irfs_opt)
        cls.obs_loc = cls.get_obs_loc(irfs_opt)
        return cls.irfs

    @staticmethod
    def get_irfs_label(required_irfs, which='both'):
        """Generate an IRF label based on options."""
        ss = 'CTAO '
        array_label = required_irfs[0].replace('SubArray', 's')
        azimuth_label = required_irfs[1].replace('AverageAz', '')
        extra = ''
        if which == 'zenith':
            extra = f' ({required_irfs[2]})'
        elif which == 'livetime':
            extra = f' ({required_irfs[3]})'
        elif which == 'both':
            extra = f' ({required_irfs[2]}-{required_irfs[3]})'
        return f"{ss}{array_label}{azimuth_label}{extra}"
    
    @staticmethod
    def get_obs_loc(irfs_opt):
        """Get observatory location based on site."""
        return observatory_locations['cta_south'] if 'South' in irfs_opt[0] else observatory_locations['cta_north']

    
    @staticmethod
    def get_irf_groups(irfs_opts):
        """Generate IRF groups with options, labels, and locations."""
        irfs_groups, irfs, irfs_labels, obs_locations = [], [], [], []
        
        arrays = irfs_opts[0] if isinstance(irfs_opts[0], list) else [irfs_opts[0]]
        azimuths = irfs_opts[1] if isinstance(irfs_opts[1], list) else [irfs_opts[1]]
        zeniths = irfs_opts[2] if isinstance(irfs_opts[2], list) else [irfs_opts[2]]
        obs_times = irfs_opts[3] if isinstance(irfs_opts[3], list) else [irfs_opts[3]]

        for array in arrays:
            for azimuth in azimuths:
                for zenith in zeniths:
                    for obs_time in obs_times:
                        irfs_opt = [array, azimuth, zenith, obs_time]
                        irfs_groups.append(irfs_opt)
                        irfs.append(Irfs.get_irfs(irfs_opt))
                        irfs_labels.append(Irfs.get_irfs_label(irfs_opt))
                        obs_locations.append(Irfs.get_obs_loc(irfs_opt))
                        
        return irfs_groups, irfs, irfs_labels, obs_locations
    
    
    @staticmethod
    def get_irfs_array(required_irfs):
        ss = 'CTAO '
        ss_0 = required_irfs[0].replace('SubArray', 's')
        ss_1 = required_irfs[1].replace('AverageAz', '')
        return  ss + ss_0 + ss_1 
    
    @staticmethod
    def get_irfs_name(required_irfs, which='both'):
        ss = 'CTAO-'
        ss_0 = required_irfs[0].replace('SubArray', 's')
        ss_1 = required_irfs[1].replace('AverageAz', '')
        _ss = ''
        if which == 'zenith':
            _ss = required_irfs[2]
        elif which == 'livetime':
            _ss = required_irfs[3]
        elif which =='both': 
            _ss = required_irfs[2] + '_' + required_irfs[3]
        return  ss + ss_0 + ss_1  + '_' + _ss
    


def calculate_annual_visibility_for_zenith_ranges(observatory_name, source_position, year=2025, time_step=30):
    """
    Calculate the annual visibility duration for a source in specific zenith angle ranges.

    Parameters:
    - observatory_name (str): The name of the observatory, as listed in `observatory_locations`.
    - source_position (SkyCoord): The celestial coordinates of the source.
    - year (int): The year for which to calculate the visibility.
    - time_step (int): Time step in minutes between each visibility calculation point.
    
    Returns:
    - annual_visibility_durations (dict): Total visibility duration in hours for each zenith angle range.
    """
    
    # Get the location of the observatory
    location = observatory_locations[observatory_name]
    
    # Define zenith angle ranges (in degrees)
    zenith_ranges = {
        '20': (10, 30),
        '40': (30, 50),
        '60': (50, 70)
    }
    # Initialize annual durations for each zenith range
    annual_visibility_durations = {range_label: 0 for range_label in zenith_ranges.keys()}
    
    # Generate time points for a typical night (18:00 to 06:00 the next morning)
    time_points = []
    for hour in range(18, 24):  # From 18:00 to 24:00
        for minute in range(0, 60, time_step):
            time_points.append(f"{hour:02d}:{minute:02d}:00")
    for hour in range(0, 6):  # From 00:00 to 06:00
        for minute in range(0, 60, time_step):
            time_points.append(f"{hour:02d}:{minute:02d}:00")
    
    # Set up the date range for the entire year
    start_date = datetime(year, 1, 1)
    end_date = datetime(year + 1, 1, 1)
    delta = timedelta(days=1)
    
    current_date = start_date
    while current_date < end_date:
        date_str = current_date.strftime("%Y-%m-%d")
        
        # Generate Time objects for all times during the night
        times = Time([f"{date_str} {time}" for time in time_points])
        
        # Create the AltAz frame for the given times and observatory location
        altaz_frame = AltAz(obstime=times, location=location)
        
        # Transform the source position to the AltAz frame to get the altitude and azimuth
        source_altaz = source_position.transform_to(altaz_frame)
        
        # Calculate zenith angles (zenith angle = 90° - altitude)
        zenith_angles = 90 * u.deg - source_altaz.alt
        
        # Count visibility duration for each zenith range
        for range_label, (zenith_min, zenith_max) in zenith_ranges.items():
            # Check if zenith angles fall within the range
            visible_times = np.sum((zenith_angles >= zenith_min * u.deg) & (zenith_angles < zenith_max * u.deg))
            
            # Convert time steps to hours and add to the annual total for this zenith range
            annual_visibility_durations[range_label] += visible_times * (time_step / 60)
        
        # Move to the next day
        current_date += delta
    
    return annual_visibility_durations
