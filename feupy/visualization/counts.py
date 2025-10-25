# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Counts class."""

import matplotlib.pyplot as plt # A collection of command style functions
from gammapy.utils.scripts import make_path
from astropy import units as u


__all__ = [
    "show_hist_counts",
]


def show_hist_counts(table, file_path=None):
    """
    Display histograms for key statistical columns in the given table.

    This function generates a figure with four subplots, each showing a histogram 
    for different statistical metrics: 
    - "counts": Number of counts detected.
    - "counts_off": Background counts.
    - "excess": Excess counts (signal-background).
    - "sqrt_ts": Significance in terms of standard deviations.

    Parameters
    ----------
    table : astropy.table.Table
        Table containing the statistical data with columns: "counts", "counts_off",
        "excess", and "sqrt_ts".
    file_path : str or Path, optional
        If provided, saves the histogram plot to the specified file path.

    Returns
    -------
    None
    """
    
    fix, axes = plt.subplots(1, 4, figsize=(12, 4))
    axes[0].hist(table["counts"])
    axes[0].set_xlabel("Counts")
    axes[0].set_ylabel("Frequency");

    axes[1].hist(table["counts_off"])
    axes[1].set_xlabel("Counts Off");

    axes[2].hist(table["excess"])
    axes[2].set_xlabel("excess");

    axes[3].hist(table["sqrt_ts"])
    axes[3].set_xlabel(r"significance ($\sigma$)");
#     path_file =  utl.get_path_counts(region_of_interest)  
#     file_name = utl.name_to_txt(file_name)
    if file_path:
        plt.savefig(file_path, bbox_inches='tight')
        
    return
