# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""ROI Map class."""

import astropy.units as u
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion, PointSkyRegion
from gammapy.maps import RegionGeom
from feupy.visualization.styles.markers import generate_catalog_markers
from gammapy.utils.scripts import make_path

__all__ = ["ROIMapPlotter"]


class ROIMapPlotter:
    """
    Class to plot a sky map with a Region of Interest (ROI) and optional sources.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        Sky coordinate of the ROI center.
    radius : `~astropy.units.Quantity`
        Radius of the ROI in angular units.
    ax : `~matplotlib.axes.Axes`, optional
        Matplotlib axes object. Default is None.
    """
    
#     def __init__(self, center, radius, ax=None):
    def __init__(self, center, radius):

        self.center = center
        self.radius = radius
#         self.ax = ax if ax is not None else plt.gca()
        self.ax = None

    def customize_legend(self, kwargs_legend):
        """
        Customize and add the legend to the plot.

        Parameters
        ----------
        kwargs_legend : dict, optional
            Dictionary defining the legend placement and style.
        """
        
        if kwargs_legend is None:
            kwargs_legend = dict(
                bbox_to_anchor=(0, -0.45),
                ncol=3,
                loc='lower left',
                markerscale=0.75,
                fontsize=5,
                labelcolor="black",
            )
        
        self.ax.legend(**kwargs_legend)
        
    def plot_roi(self, color='blue', linestyle='--'):
        """Plot the Region of Interest (ROI) on the sky map."""
        region = RegionGeom(CircleSkyRegion(self.center, self.radius))
        self.ax = region.plot_region(color=color, linestyle=linestyle)
        return self.ax

    def plot_sources(self, sources, ref_markers=None):
        """
        Plot astronomical sources on the sky map.

        Parameters
        ----------
        sources : list
            List of sources with position and name.
        ref_markers : dict, optional
            Dictionary of reference markers for each source.
        """
        ref_markers = ref_markers or generate_catalog_markers(sources, marker_size=6, PALETTE=None)
            
        # Plot each source with corresponding marker
        for index, source in enumerate(sources):
            
            kwargs_point = ref_markers[sources.labels[index]]     
            kwargs_point.update({'fillstyle': 'full', 'lw': 0})

            point = RegionGeom(PointSkyRegion(center=source.position))
            self.ax = point.plot_region(
                ax=self.ax,
                facecolor=kwargs_point['color'],
                edgecolor='black',
                kwargs_point=kwargs_point,
            )
            


    def set_axes(self, xlabel="R.A. (J2000)", ylabel="Dec. (J2000)", size=12):
        """
        Set the axis labels and font size for the plot.

        Parameters
        ----------
        xlabel : str, optional
            Label for the x-axis. Default is 'R.A. (J2000)'.
        ylabel : str, optional
            Label for the y-axis. Default is 'Dec. (J2000)'.
        size : int, optional
            Font size for the labels. Default is 12.
        """
        self.ax.set_xlabel(xlabel, size=size)
        self.ax.set_ylabel(ylabel, size=size)
        self.ax.grid(True)

    def add_roi_text(self):
        """Add text annotation for the ROI on the plot."""
        self.ax.text(0.1, 0.93, f'ROI ({self.radius})', transform=self.ax.transAxes)

    def save_plot(self, file_path):
        """
        Save the sky map to a file.

        Parameters
        ----------
        file_path : str or `~pathlib.Path`
            Path to save the plot.
        """
        if file_path:
            plt.savefig(file_path, bbox_inches='tight')

    def plot(self, sources=None, file_path=None, **kwargs):
        """
        Display the ROI map with optional sources and save it to a file.

        Parameters
        ----------
        sources : list of sources, optional
            List of sources with position and name.
        file_path : str or `~pathlib.Path`, optional
            Path to save the plot. Default is None.
        **kwargs : dict, optional
            Additional keyword arguments for customization.
        """
        
        # Use the provided axes or get the current axes
        self.ax = plt.gca() if  self.ax is None else  self.ax
        
        # Plot ROI
        self.plot_roi()

        # Plot sources if provided
        if sources:
            ref_markers = kwargs.get('ref_markers')
            self.plot_sources(sources, ref_markers)

        # Customize legend
        kwargs_legend = kwargs.get('kwargs_legend', None)
        self.customize_legend(kwargs_legend)

        # Set axis labels
        self.set_axes()

        # Add ROI text
        self.add_roi_text()

        # Save the map if file path is provided
        self.save_plot(file_path)

        return self.ax
    
    
# def add_annotations(self, annotations, **kwargs):
#     """Add annotations to the ROI map."""
#     # Function implementation...

#     kwargs.setdefault('ref_markers', {})
#     ref_markers = kwargs['ref_markers']


# def customize_display(self, styles=None, **kwargs):
#     """Customize the appearance of the ROI map."""
#     # Function implementation...

# def save_map(self, file_path):
#     """Save the map to a file."""
#     if self.ax:
#         self.ax.figure.savefig(file_path, bbox_inches='tight')