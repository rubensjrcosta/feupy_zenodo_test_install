# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Sensitivity class."""

from gammapy.maps.axes import UNIT_STRING_FORMAT
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt
import astropy.units as u
from feupy.analysis.irfs import Irfs
from feupy.visualization.utils.units import DEFAULT_UNIT, DEFAULT_XAXIS_LABEL, DEFAULT_YAXIS_LABEL


__all__ = ["plot_irfs"]


def plot_irfs(
    tables_south,
    model=None, energy_bounds=None, ax=None, which='e2dnde', which_label='both', 
    int_sens_label=True, fname=None
):
    """
    Plot IRFs for both South and North tables, with optional model overlay.

    Parameters
    ----------
    tables_south : list of Table
        List of tables for the southern observations.
    tables_north : list of Table
        List of tables for the northern observations.
    model : SpectralModel, optional
        Spectral model to overlay on the plot.
    energy_bounds : tuple, optional
        Energy bounds for model plotting.
    ax : matplotlib.axes.Axes, optional
        Axis on which to plot. A new axis is created if not provided.
    which : str, default='e2dnde'
        Y-axis variable for sensitivity plot.
    which_label : str, default='both'
        Type of IRF label to display.
    int_sens_label : bool, default=True
        Whether to display integrated sensitivity label.
    fname : str, optional
        File name to save the plot, if provided.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The plot axis.
    """
    if not isinstance(tables_south, list) :
        raise TypeError("Both `tables_south` and `tables_north` must be lists of tables.")
    
    linestyle = ['solid', 'solid', 'solid']
    #linestyle = ['solid', (0, (5, 1)), (0, (3, 5, 1, 5))]
    ax = plt.gca() if ax is None else ax

    ax.set_prop_cycle(color=['blue', 'red', 'green'], linestyle=linestyle)
    for index, table in enumerate(tables_south):
        required_irfs = table.meta['IRFS']
        label = Irfs.get_irfs_label(required_irfs, which=which_label)
        if int_sens_label:
            int_sens = u.Quantity(table.meta['INT_SENS'])
            unit = int_sens.unit.to_string(UNIT_STRING_FORMAT)
            label += f' ({int_sens.value:.2e} {unit})'
        ax = plot_sensitivity_from_table(table, which=which, ax=ax, label=label)

#     ax.set_prop_cycle(color=['green', 'green', 'green'], linestyle=linestyle)
#     for index, table in enumerate(tables_north):
#         required_irfs = table.meta['IRFS']
#         label = Irfs.get_irfs_label(required_irfs, which=which_label)
#         if int_sens_label:
#             int_sens = u.Quantity(table.meta['INT_SENS'])
#             unit = int_sens.unit.to_string(UNIT_STRING_FORMAT)
#             label += f' ({int_sens.value:.2e} {unit})'
#         ax = plot_sensitivity_from_table(table, which=which, ax=ax, label=label)

    if model is not None:
        kwargs = {'ax': ax, 'sed_type': 'e2dnde'}
        ax = model.spectral_model.plot(energy_bounds=energy_bounds, label=model.name, color="k", **kwargs)
        ax = model.spectral_model.plot_error(energy_bounds=energy_bounds, **kwargs)

    ax.legend(loc="best", scatterpoints=1, handlelength=3, fontsize=6)

    return ax


# def plot_tables_sensitivity(
#     tables_south, tables_north, model=None, sens_info=None, box_name=None,
#     abs_model=None, which_label='both', file_path=None, **kwargs
# ):
#     """Plot sensitivity tables with optional models and save to file if specified."""

#     if not isinstance(tables_south, list) or not isinstance(tables_north, list):
#         raise TypeError("Both `tables_south` and `tables_north` must be lists of tables.")
        
#     kwargs.setdefault('energy_bounds', [3e-2, 1e2] * u.TeV)
#     energy_bounds = kwargs['energy_bounds']

#     kwargs.setdefault('ylim', [1e-14, 1e-8])
#     ylim = kwargs['ylim']

#     kwargs.setdefault('sed_type', 'e2dnde')

#     sed_type = kwargs['sed_type']

#     fig, ax = plt.subplots()

#     ax = plot_irfs_superpose(
#         tables_south, tables_north, ax=ax, which='e2dnde', 
#         which_label=which_label, int_sens_label=False
#     )

#     if model is not None:
#         model.spectral_model.plot(label =  model.name,
#                                   energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), 
#                                   ax=ax, 
#                                   sed_type = sed_type,
#                                   linestyle='solid', marker=',', color='black')
#         model.spectral_model.plot_error(energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), sed_type = sed_type, ax=ax)

#         if abs_model is not None:
#             abs_model.spectral_model.plot(label =  abs_model.name,
#                                           energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), ax=ax, 
#                                                                             sed_type = kwargs['sed_type'],
#                                           linestyle='--', marker=',', color='black')
#             abs_model.spectral_model.plot_error(energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), sed_type = kwargs['sed_type'],ax=ax)

#     ax.set_xlabel(kwargs.get('xaxis_label', DEFAULT_XAXIS_LABEL['TeV']))
#     ax.set_ylabel(kwargs.get('yaxis_label', DEFAULT_YAXIS_LABEL[kwargs['sed_type']]))
#     ax.set_ylim(ylim)
#     ax.set_xlim(energy_bounds.value)
#     ax.tick_params(which="both", labelbottom=True, labeltop=False, labelleft=True, labelright=False,
#                  bottom=True, top=True, left=True, right=True,direction="in")
    
#     ax.legend(loc='upper right')
#     if box_name:
#         ax.text(0.1, 0.9, box_name, transform=ax.transAxes)
#     if sens_info:
#         ax.text(.1, .05, sens_info, fontsize=6, transform=ax.transAxes)
    
#     if file_path:
#         plt.savefig(file_path)
#         print(f"Sensitivity plot saved to {file_path}.")

#     return fig, ax

def plot_sensitivity_from_table(
    sens_table, which='e2dnde', ax=None, plot_xerr=False, **kwargs
):
    """ 
    Plot sensitivity or related metrics from a table on a given axis.

    Parameters
    ----------
    sens_table : Table
        Table containing sensitivity data. Expected columns include 'e_ref', 
        and one of 'e2dnde', 'excess', 'background', or 'on_radii' depending on `which`.
    which : str, default='e2dnde'
        Specifies which metric to plot. Options are 'e2dnde' (flux sensitivity), 
        'excess' (excess counts), 'background' (background counts), or 'on_radii' (on-region radius).
    ax : matplotlib.axes.Axes, optional
        Axis on which to plot. A new axis is created if not provided.
    plot_xerr : bool, default=False
        Whether to plot horizontal error bars based on 'e_min' and 'e_max' columns.
    **kwargs : dict
        Additional keyword arguments passed to `ax.errorbar` for plot customization.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axis with the plotted sensitivity.
    """
    ax = plt.gca() if ax is None else ax
    e = sens_table["e_ref"]
    s = sens_table[which]

    xlabel = f"Energy [{e.unit.to_string(UNIT_STRING_FORMAT)}]"
    ylabel = {
        'excess': 'Excess counts',
        'background': "Background counts",
        'on_radii': f"On region radius [{s.unit.to_string(UNIT_STRING_FORMAT)}]",
        'e2dnde': f"Flux Sensitivity [{s.unit.to_string(UNIT_STRING_FORMAT)}]"
    }.get(which, f"{which} [{s.unit.to_string(UNIT_STRING_FORMAT)}]")

    xerr = (sens_table["e_max"] - sens_table["e_min"]) / 2 if plot_xerr else None

    with quantity_support():
        ax.errorbar(e, s, xerr=xerr, **kwargs)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="lower left", scatterpoints=1, handlelength=3, fontsize=8)

    return ax


def plot_irfs_superpose(
    tables_south, tables_north,
    model=None, energy_bounds=None, ax=None, which='e2dnde', which_label='both', 
    int_sens_label=True, fname=None
):
    """
    Plot IRFs for both South and North tables, with optional model overlay.

    Parameters
    ----------
    tables_south : list of Table
        List of tables for the southern observations.
    tables_north : list of Table
        List of tables for the northern observations.
    model : SpectralModel, optional
        Spectral model to overlay on the plot.
    energy_bounds : tuple, optional
        Energy bounds for model plotting.
    ax : matplotlib.axes.Axes, optional
        Axis on which to plot. A new axis is created if not provided.
    which : str, default='e2dnde'
        Y-axis variable for sensitivity plot.
    which_label : str, default='both'
        Type of IRF label to display.
    int_sens_label : bool, default=True
        Whether to display integrated sensitivity label.
    fname : str, optional
        File name to save the plot, if provided.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The plot axis.
    """
    if not isinstance(tables_south, list) or not isinstance(tables_north, list):
        raise TypeError("Both `tables_south` and `tables_north` must be lists of tables.")
    
    linestyle = ['solid', (0, (5, 1)), (0, (3, 5, 1, 5))]
    ax = plt.gca() if ax is None else ax

    ax.set_prop_cycle(color=['blue', 'blue', 'blue'], linestyle=linestyle)
    for index, table in enumerate(tables_south):
        required_irfs = table.meta['IRFS']
        label = Irfs.get_irfs_label(required_irfs, which=which_label)
        if int_sens_label:
            int_sens = u.Quantity(table.meta['INT_SENS'])
            unit = int_sens.unit.to_string(UNIT_STRING_FORMAT)
            label += f' ({int_sens.value:.2e} {unit})'
        ax = plot_sensitivity_from_table(table, which=which, ax=ax, label=label)

    ax.set_prop_cycle(color=['green', 'green', 'green'], linestyle=linestyle)
    for index, table in enumerate(tables_north):
        required_irfs = table.meta['IRFS']
        label = Irfs.get_irfs_label(required_irfs, which=which_label)
        if int_sens_label:
            int_sens = u.Quantity(table.meta['INT_SENS'])
            unit = int_sens.unit.to_string(UNIT_STRING_FORMAT)
            label += f' ({int_sens.value:.2e} {unit})'
        ax = plot_sensitivity_from_table(table, which=which, ax=ax, label=label)

    if model is not None:
        kwargs = {'ax': ax, 'sed_type': 'e2dnde'}
        ax = model.spectral_model.plot(energy_bounds=energy_bounds, label=model.name, color="k", **kwargs)
        ax = model.spectral_model.plot_error(energy_bounds=energy_bounds, **kwargs)

    ax.legend(loc="best", scatterpoints=1, handlelength=3, fontsize=6)

    return ax


def plot_tables_sensitivity(
    tables_south, tables_north=None, model=None, sens_info=None, box_name=None,
    abs_model=None, which_label='both', file_path=None, **kwargs
):
    """Plot sensitivity tables with optional models and save to file if specified."""

    if not isinstance(tables_south, list):
        raise TypeError("`tables` must be lists of tables.")
    if tables_north:
        if not isinstance(tables_north, list):
            raise TypeError("`tables` must be lists of tables.")

        
    kwargs.setdefault('energy_bounds', [3e-2, 1e2] * u.TeV)
    energy_bounds = kwargs['energy_bounds']

    kwargs.setdefault('ylim', [1e-14, 1e-8])
    ylim = kwargs['ylim']

    kwargs.setdefault('sed_type', 'e2dnde')

    sed_type = kwargs['sed_type']

    fig, ax = plt.subplots()
    if tables_north:

        ax = plot_irfs_superpose(
            tables_south, tables_north, ax=ax, which='e2dnde', 
            which_label=which_label, int_sens_label=False
        )
    else: 
        ax = plot_irfs(
            tables_south, ax=ax, which='e2dnde', 
            which_label=which_label, int_sens_label=False
        )

    if model is not None:
        model.spectral_model.plot(label =  model.name,
                                  energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), 
                                  ax=ax, 
                                  sed_type = sed_type,
                                  linestyle='solid', marker=',', color='black')
        model.spectral_model.plot_error(energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), sed_type = sed_type, ax=ax)

        if abs_model is not None:
            abs_model.spectral_model.plot(label =  abs_model.name,
                                          energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), ax=ax, 
                                                                            sed_type = kwargs['sed_type'],
                                          linestyle='--', marker=',', color='black')
            abs_model.spectral_model.plot_error(energy_bounds=kwargs.get('energy_bounds_fit', energy_bounds), sed_type = kwargs['sed_type'],ax=ax)

    ax.set_xlabel(kwargs.get('xaxis_label', DEFAULT_XAXIS_LABEL['TeV']))
    ax.set_ylabel(kwargs.get('yaxis_label', DEFAULT_YAXIS_LABEL[kwargs['sed_type']]))
    ax.set_ylim(ylim)
    ax.set_xlim(energy_bounds.value)
    ax.tick_params(which="both", labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                 bottom=True, top=True, left=True, right=True,direction="in")
    
    ax.legend(loc='upper right', frameon=False)
    if box_name:
        ax.text(0.1, 0.9, box_name, transform=ax.transAxes)
    if sens_info:
        ax.text(.1, .05, sens_info, fontsize=6, transform=ax.transAxes)
    
    if file_path:
        plt.savefig(file_path)
        print(f"Sensitivity plot saved to {file_path}.")

    return fig, ax
