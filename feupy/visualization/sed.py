# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
SEDPlotter class.

This module provides the SEDPlotter class, which is used to plot Spectral Energy Distributions (SEDs)
from a collection of datasets and models. It offers flexible options for customizing plot appearance,
legends, axis labels, units, and plot limits.
"""

import matplotlib.pyplot as plt
from astropy import units as u
from feupy.visualization.styles.markers import generate_specified_marker_set, get_linestyles
from feupy.utils.datasets import get_energy_bounds_from_datasets
from feupy.visualization.utils.units import DEFAULT_UNIT, DEFAULT_XAXIS_LABEL, DEFAULT_YAXIS_LABEL 
from feupy.visualization import LINESTYLES_DEFAULT


__all__ = ["SEDPlotter"]


class SEDPlotter:
    def __init__(self, datasets, models=None, sed_type="e2dnde"):
        """
        Initialize the SEDPlotter with datasets and models.

        Parameters
        ----------
        datasets : `~gammapy.datasets.Datasets`
            Collection of datasets from which to plot SED data.
        models : `~gammapy.modeling.models.Models`, optional
            Models that describe the spectral energy distribution (SED). Default is None.
        sed_type : {"dnde", "e2dnde"}, optional
            Type of SED plot to generate. Default is "e2dnde".
        """
        self.datasets = datasets
        self.models = models
        self.sed_type = sed_type
        self.ax = None
        
    def customize_legend(self, **kwargs):
        """
        Configure the legend on the plot with customizable placement.

        Parameters
        ----------
        kwargs : dict
            Dictionary for legend placement settings.
        """
        kwargs.setdefault(
            'kwargs_legend',
            dict(
                # bbox_to_anchor=(0, -0.45),
                ncol=3,
                loc='lower left',
                markerscale=0.75,
                fontsize=5,
                labelcolor="black",
                frameon=False,
            )
        )
        self.ax.legend(**kwargs['kwargs_legend'])
        
    def set_axis_labels(self, **kwargs):
        """
        Set axis labels and units for the plot.

        Parameters
        ----------
        kwargs : dict
            Dictionary containing settings for axis labels and units.
        """

        xlabel, ylabel = kwargs['axis']['label']        
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        
    def set_axis_units(self, **kwargs):
        """
        Set axis labels and units for the plot.

        Parameters
        ----------
        kwargs : dict
            Dictionary containing settings for axis labels and units.
        """
        xunits, yunits = kwargs['axis']['units']
        self.ax.xaxis.set_units(u.Unit(xunits))
        self.ax.yaxis.set_units(u.Unit(yunits))
        
    def set_plot_limits(self, **kwargs):
        """
        Set the plot limits for the energy and flux axes.

        Parameters
        ----------
        kwargs : dict
            Dictionary containing settings for plot limits.
        """
        energy_bounds = kwargs['limits']['energy_bounds']
        ylim = kwargs['limits']['ylim']
        self.ax.set_xlim(energy_bounds.value)
        self.ax.set_ylim(ylim)

    def plot_datasets(self, ax, ref_markers, plot_kwargs):
        """
        Plot the datasets on the given axis.

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            Matplotlib axes object.
        ref_markers : dict
            Marker settings for the datasets.
        plot_kwargs : dict
            Additional keyword arguments for plotting.
        """
        for dataset in self.datasets:
            kwargs_dataset = ref_markers.get(dataset.name, {}).copy()

            kwargs_dataset.update(
                dict(ls='None', lw=0.5, markeredgecolor='k', mew=0.8, elinewidth=0.6, capsize=1.5)
            )
            
            # Plot dataset data
            dataset.data.plot(**plot_kwargs, **kwargs_dataset)

            # Prepare model plot arguments
            color = kwargs_dataset.pop('color', 'black')
            energy_bounds = get_energy_bounds_from_datasets(dataset)
            kwargs_model = dict(edgecolor=color, alpha=0.2, color='black', facecolor=color, 
                                energy_bounds=energy_bounds)
            
            # Plot model error
            if dataset.models is not None and dataset.name in dataset.models.names:

#             dataset.models[0].spectral_model.plot_error(**plot_kwargs, **kwargs_model)
                dataset.models[dataset.name].spectral_model.plot_error(**plot_kwargs, **kwargs_model)

    def plot_models(self, ax, ref_markers, plot_kwargs, energy_bounds, error_band= False):
        """
        Plot the models if provided.

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            Matplotlib axes object.
        kwargs_models : dict
            Additional keyword arguments for plotting models.
        plot_kwargs : dict
            Additional keyword arguments for plotting.
        energy_bounds : `~astropy.units.Quantity`
            Energy bounds for the model plot.
        """
        if self.models:            
            linestyles = get_linestyles(LINESTYLES_DEFAULT, num_linestyles=len(self.models))
            for index, model in enumerate(self.models):

                # Prepare model plot arguments
                kwargs_model = dict(
                    label=model.name, 
                    marker=',', 
                    color='black',
                    linestyle=linestyles[index], 
                    energy_bounds=energy_bounds,
                )
    
                spectral_model = model.spectral_model
                spectral_model.plot(**plot_kwargs, **kwargs_model)
                if error_band:
                    spectral_model.plot_error(energy_bounds=energy_bounds,alpha=0.05,**plot_kwargs)

    def save_plot(self, file_path):
        """
        Save the SED to a file.

        Parameters
        ----------
        file_path : str or `~pathlib.Path`
            File path or name where the plot will be saved.
        """
        if file_path:
            plt.savefig(file_path, bbox_inches='tight')

    def plot(self, ax=None, file_path=None, ref_markers=None,  box_name=None, error_band=False, **kwargs):
        """
        Generate and display the SED plot.

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`, optional
            Matplotlib axes object to draw the plot on.
        file_path : str or `~pathlib.Path`, optional
            File path or name where the plot will be saved.
        ref_markers : dict, optional
            Marker settings for the datasets.
        **kwargs : dict, optional
            Additional settings for plot customization.
        """
        # Generate marker set for datasets
        refs_names = self.datasets.names
        if self.models:
            refs_names.extend(self.models.names)
        ref_markers = ref_markers or generate_specified_marker_set(refs_names, marker_size=4)

        # Use the provided axes or get the current axes
        self.ax = ax if ax else plt.gca()

        # Set axis units and label
        kwargs.setdefault(
            'axis',
            dict(
                label=(DEFAULT_XAXIS_LABEL['TeV'], DEFAULT_YAXIS_LABEL[self.sed_type]),
                units=('TeV', 'TeV cm-2 s-1'),
            )
        )
            
        self.set_axis_units(**kwargs)

        # Plot configuration settings
        plot_kwargs = {
            "ax": self.ax,
            "sed_type": self.sed_type,
        }


        # Plot the models if provided
        kwargs.setdefault('kwargs_models', {})
        kwargs_models = kwargs['kwargs_models']
    
        kwargs.setdefault(
            'limits',
            dict(
                energy_bounds=[1e-5, 2e3] * u.TeV,
                ylim=[1e-23, 1e-7]
            )
        )
#         print(kwargs_models)
        # Plot the datasets
        self.plot_datasets(self.ax, ref_markers, plot_kwargs)
        
        if kwargs_models:
            if 'energy_bounds' in list(kwargs_models.keys()):
                energy_bounds = kwargs_models['energy_bounds']
        else: energy_bounds = kwargs['limits']['energy_bounds']
        self.plot_models(self.ax, ref_markers, plot_kwargs, energy_bounds, error_band)
    
        # Set plot limits
        self.set_plot_limits(**kwargs)

        # Configure legend
        self.customize_legend(**kwargs)
#         ax.legend(**leg_place)

        if box_name:
            self.ax.text(0.1, 0.9, box_name, transform=self.ax.transAxes)
        
        # Set axis labels
        self.set_axis_labels(**kwargs)
        
        # Save plot if file_path is provided
        self.save_plot(file_path)

        return self.ax 
    
    
# from astropy import units as u
# import matplotlib.pyplot as plt # A collection of command style functions

# def SED_from_catalogs(
#     counterparts, datasets_counterparts, models_counterparts, colors_dict, region_of_interest,
#                            sed_type = "e2dnde", 
#                            axis_dict =  dict(
#     label =  (r'$\rm{E\ [TeV] }$', r'$\rm{E^2\ J(E)\ [TeV\ cm^{-2}\ s^{-1}] }$'),
#     units =  (          'TeV',                       'TeV  cm-2     s-1')
# ),                         
#     energy_bounds = [1e-5, 1e2] * u.TeV, 
#     ylim = [1e-13, 1e-9]
#                           ):
    
#     sed_type = sed_type
#     axis_dict =axis_dict
#     energy_bounds =energy_bounds
#     ylim = ylim
    
        
#     for index, (counterpart, dataset, model) in enumerate(zip(counterparts, datasets_counterparts, models_counterparts)):    

#         counterpart_name = counterpart.name
#         flux_points = counterpart.flux_points
# #         spectral_model = model.spectral_model
#         spectral_model = counterpart.spectral_model()
#         spectral_model_tag = spectral_model.tag[0]
#         spectral_model_tag_short = spectral_model.tag[1]

#         ax = plt.subplot()
#         ax.xaxis.set_units(u.Unit(axis_dict['units'][0]))
#         ax.yaxis.set_units(u.Unit(axis_dict['units'][1]))
    
#         xlabel = axis_dict['label'][0]
#         ylabel = axis_dict['label'][1]

#         kwargs = {
#             "ax": ax, 
#             "sed_type": sed_type
#         }
#         kwargs_fit = {
#             "label": f"{spectral_model_tag_short} (fit)"
#         }

#         color = colors_dict[dataset.name][0]
#         marker = colors_dict[dataset.name][1]
#         flux_points.plot(label = dataset.name, color= color, marker = marker, **kwargs)

#         energy_bounds = flux_points.energy_min[0], flux_points.energy_max[-1]
#     #     e2dnde_errn = flux_points.e2dnde_errn.data
#     #     e2dnde_errp = flux_points.e2dnde_errp.data
#     #     ylim =min(np.nanmin(e2dnde_errp),np.nanmin(e2dnde_errn)), max(np.nanmax(e2dnde_errp),np.nanmax(e2dnde_errn))

#         spectral_model.plot(energy_bounds=energy_bounds, marker = ',', ls= "-", color="k", **kwargs, **kwargs_fit)
#         spectral_model.plot_error(energy_bounds=energy_bounds, **kwargs)

#     #     kwargs_spectrum = {"kwargs_model": {"color":"red", "ls":"--"}, "kwargs_fp":{"color":"green", "marker":"o"}}
#     #     dataset_fp.plot_spectrum(**kwargs, **kwargs_spectrum)  
#     #    ax.set_ylim(ylim)
    
#         ax.set_xlim(energy_bounds)
#         plt.xlabel(xlabel)   
#         plt.ylabel(ylabel)
#         plt.legend()
#         file_name = f'{utl.name_to_txt(dataset.name.replace(":",""))}_{spectral_model_tag_short}{cfg.format_png}'
#         file_path = utl.get_path_SED_from_catalogs(region_of_interest) / file_name 
#         plt.savefig(file_path, bbox_inches='tight')
#         plt.show()
