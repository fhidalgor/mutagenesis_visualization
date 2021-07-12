"""
This module contains the class to plot a regular heatmap from enrichment
scores.
"""
from typing import Any, Dict, Union
from pathlib import Path
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

from mutagenesis_visualization.main.default_parameters.kwargs import default_kwargs
from mutagenesis_visualization.main.default_parameters.font_parameters import font_parameters
from mutagenesis_visualization.main.heatmaps.heatmap_utils import (
    generate_cartoon,
    labels,
    hierarchical_sort,
)
from mutagenesis_visualization.main.utils.snv import add_snv_boolean
from mutagenesis_visualization.main.utils.pandas_functions import df_rearrange
from mutagenesis_visualization.main.utils.save_work import save_work


def plot_heatmap(
    self,
    nancolor: str = 'lime',
    show_cartoon: bool = False,
    show_snv: bool = False,
    hierarchical: bool = False,
    output_file: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
) -> None:
    '''
    Generate a heatmap plot of the enrichment scores.

    Parameters
    ----------
    self : object from class *Screen*

    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.

    show_carton : boolean, default False
        If true, the plot will display a cartoon with the secondary structure.
        The user must have added the secondary structure to the object.

    show_snv : boolean, default False
        If true, it will only display mutants that are a single nucleotide
        variant (SNV) of the wild-type protein sequence. The algorithm
        does not take into account the wild-type DNA allele, so it will
        include any possible mutant that is one base away.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax_object).
    Returns
    ----------
    fig, ax_object, ax_object2, ax_object3, averageresidue : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    '''

    # load font parameters
    font_parameters()

    # update kwargs
    temp_kwargs: Dict[str, Any] = copy.deepcopy(default_kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]

    # sort data by rows in specified order by user
    df_data: pd.DataFrame = df_rearrange(
        add_snv_boolean(self.dataframe_stopcodons.copy()),
        temp_kwargs['neworder_aminoacids'],
        values='Score_NaN',
        show_snv=show_snv
    )

    # average of residues by positon
    average = [add_snv_boolean(self.dataframe.copy()).groupby(by='Position').mean()['Score_NaN']]
    # For 1 column case

    # Create new sequence that we may to change order later
    self.sequence_updated = self.sequence

    # sort data by columns if specified by user
    if hierarchical:
        sorted_columns = hierarchical_sort(df_data)
        # adds the initial number so index match
        sorted_columns_corrected = sorted_columns + self.start_position
        df_data = df_data[sorted_columns_corrected]
        temp = pd.DataFrame(average)  # had to convert back to df_data to be able to sort properly
        average = temp[sorted_columns_corrected]
        # Change the sequence order
        self.sequence_updated = pd.DataFrame(list(self.sequence)).T[sorted_columns]
        self.sequence_updated = list(self.sequence_updated.iloc[0])
    elif len(average) == 1:
        average = [np.array(average[0])]

    # declare figure and subplots
    figwidth = 14 * len(df_data.columns) / 165
    temp_kwargs['figsize_x'] = kwargs.get('figsize_x', figwidth)
    # Change parameters depending on whether cartoon is on or off
    if show_cartoon:
        figheight = 2.45
        temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)
        fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
        gs_object = gridspec.GridSpec(
            nrows=3,
            ncols=2,
            height_ratios=[len(df_data), 1, 5],
            width_ratios=[len(df_data.columns), 1],
        )
    else:
        figheight = 2
        temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)
        fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
        gs_object = gridspec.GridSpec(
            nrows=2,
            ncols=2,
            height_ratios=[len(df_data), 1],
            width_ratios=[len(df_data.columns), 1],
        )

    ax_object = plt.subplot(gs_object[0, 0])
    averageresidue = plt.subplot(gs_object[1, 0])
    cbar1 = plt.subplot(gs_object[0, 1])
    cbar2 = plt.subplot(gs_object[1, 1])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax_object.pcolormesh(
        df_data,
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=cmap,
        edgecolors='k',
        linewidths=0.2,
        antialiased=True,
        color='darkgrey'
    )

    # average by position
    averageresidue.pcolormesh(
        average,
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=cmap,
        edgecolors='k',
        linewidths=0.2,
        antialiased=True,
        color='darkgrey'
    )

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax_object.set_xticks(np.arange(len(df_data.columns)) + 0.5, minor=False)
    ax_object.set_yticks(np.arange(len(df_data)) + 0.5, minor=False)

    # position of axis labels
    ax_object.tick_params('x', direction='out', pad=-2.5)
    ax_object.tick_params('y', direction='out', pad=0.4)

    # make new axes
    ax_object2 = ax_object.twiny()
    ax_object3 = ax_object.twinx()

    # tune the axes
    ax_object2.set_xticks(np.arange(len(df_data.columns)) + 0.5, minor=False)
    ax_object3.set_yticks(np.arange(len(df_data)) + 0.5, minor=False)
    ax_object2.tick_params(direction='out', pad=4)
    ax_object3.tick_params(direction='out', pad=0.4)
    averageresidue.tick_params(direction='out', pad=-2)
    averageresidue.set_xticks(
        np.arange(len(df_data.columns)) + 0.5,
        minor=False,
    )
    averageresidue.set_yticks(np.arange(0.5) + 0.5)

    # Set the limits of the new axis from the original axis limits
    ax_object2.set_xlim(ax_object.get_xlim())
    ax_object3.set_ylim(ax_object.get_ylim())

    # want a more natural, table-like display
    ax_object.invert_yaxis()
    ax_object.xaxis.tick_top()
    ax_object3.invert_yaxis()

    # remove ticks
    ax_object.xaxis.set_ticks_position('none')
    ax_object.yaxis.set_ticks_position('none')
    ax_object2.yaxis.set_ticks_position('none')
    ax_object2.xaxis.set_ticks_position('none')
    ax_object3.yaxis.set_ticks_position('none')
    averageresidue.xaxis.set_ticks_position('none')
    averageresidue.yaxis.set_ticks_position('none')

    # so labels of x and y do not show up and my labels show up instead
    ax_object.set_xticklabels(
        list(self.sequence_updated),
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False,
    )
    ax_object.set_yticklabels(
        temp_kwargs['neworder_aminoacids'],
        fontsize=6,
        fontname="Arial",
        color='k',
        minor=False,
    )
    # For numbering labels, change if hierarchical sorting is true
    if not hierarchical:
        ax_object2.set_xticklabels(
            temp_kwargs['number_sequencelabels'][0 : len(df_data.columns)],
            fontsize=10,
            fontname="Arial",
            color='k',
            minor=False
        )
    else:
        ax_object2.tick_params(direction='out', pad=7)
        ax_object2.set_xticklabels(
            sorted_columns_corrected,
            fontsize=5,
            fontname="Arial",
            color='k',
            minor=False,
            rotation=90,
        )
    ax_object3.set_yticklabels(
        temp_kwargs['neworder_aminoacids'],
        fontsize=6,
        fontname="Arial",
        color='k',
        minor=False,
    )
    averageresidue.set_xticklabels(
        list(self.sequence_updated),
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False,
    )
    rowaverage = ''
    averageresidue.set_yticklabels(rowaverage, fontsize=6, fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax_object.get_yticklabels():
        ylabel.set_horizontalalignment('center')
    for ylabel in ax_object3.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # for coloring the residues that are 10,20...
    for xtick, color in zip(ax_object.get_xticklabels(), temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)
    for xtick, color in zip(averageresidue.get_xticklabels(), temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)
    # _____________________________________________________________________________

    # for color bar format
    cbar1.axis('off')
    cbar2.axis('off')
    cb_object = plt.colorbar(
        heatmap,
        fraction=1,
        pad=0,
        ax_object=[cbar1],
        aspect=5,
        ticks=[
            temp_kwargs['colorbar_scale'][0],
            np.mean(temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]
        ],
        orientation='vertical'
    )
    cb_object.ax_object.set_yticklabels(
        cb_object.ax_object.get_yticklabels(),
        fontsize=6,
        fontname="Arial",
        color='k',
    )
    cb_object.update_ticks()
    plt.text(
        1.2 + 10 / len(df_data.columns),
        0.7,
        r'$\langle∆E^x_i\rangle_x$',
        transform=cbar1.transAxes,
        horizontalalignment='center',
        fontsize=7,
        fontname="Arial",
        color='k'
    )

    gs_object.update(hspace=0.1, wspace=0.1 / len(df_data.columns) * 50)

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=12)

    # Cartoon
    if show_cartoon:
        generate_cartoon(self, gs_object, 2, temp_kwargs['cartoon_colors'], 0.025)

    # save file
    save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax_object, ax_object2, ax_object3, averageresidue

    if temp_kwargs['show']:
        plt.show()
