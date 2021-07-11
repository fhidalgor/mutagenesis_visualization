"""

"""
import numpy as np
import pandas as pd
import itertools
import copy
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from pathlib import Path
from typing import Union

def plot_position(
    self, position, output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Choose a position and plot in a bargraph the enrichment score for each
    substitution. Red for gain of function, blue for loss of function.

    Parameters
    ----------
    self : object from class *Screen*

    position : int
        number of residue of the protein to display.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).

    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # load parameters
    _parameters_mean()

    # Select position
    df = self.dataframe.loc[self.dataframe['Position'] == position].copy()

    # Color
    df['Color'] = df.apply(
        code_utils._color_data,
        axis=1,
        args=(temp_kwargs['color_gof'], temp_kwargs['color_lof'])
    )

    # make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    width = 0.5

    # Color based on values
    ax.bar(df['Aminoacid'], df['Score'], width, color=df['Color'], ec='k')

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(
        temp_kwargs['y_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=10,
        rotation=0
    )

    ax.set_xlabel(
        'Residue', fontsize=10, fontname="Arial", color='k', labelpad=4
    )
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()
