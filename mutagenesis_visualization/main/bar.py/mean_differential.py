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

def plot_meandifferential(
    self,
    obj2,
    show_cartoon=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    Plot the mean positional difference between two experiments.

    Parameters
    ----------
    self : object from class *Screen*

    obj2 : another Screen object to compare with

    show_cartoon : boolean, default False
        If true, the plot will display a cartoon with the secondary structure. The user must have added the secondary structure to the object.

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
    temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2.5))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
    temp_kwargs['y_label'] = kwargs.get(
        'y_label', r'Mean Differential $∆E^i_x$'
    )

    # load parameters
    _parameters_mean()

    # make pandas
    df = code_utils._process_meanresidue(self, obj2)

    # make cartoon
    if show_cartoon:
        fig = plt.figure(figsize=temp_kwargs['figsize'])
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
        ax = plt.subplot(gs[0])
    else:
        fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # plot
    ax.plot(df['Position'], df['d1 - d2'], color='k')

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(
        temp_kwargs['y_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=-5,
        rotation=90
    )
    ax.set_xticks(
        np.arange(self.start_position,
                  len(df) + self.start_position, 20)
    )
    ax.set_xlabel(
        'Residue', fontsize=10, fontname="Arial", color='k', labelpad=4
    )
    ax.set_xlim(
        self.start_position - 0.1,
        len(df) + self.start_position - 1 + 0.1
    )

    # cartoon
    title_pad = 0
    if show_cartoon:
        title_pad = 2.5
        obj = obj2
        if len(self.dataframe) < len(obj2.dataframe):
            obj = self
            _generate_cartoon(
                obj,
                gs,
                1,
                temp_kwargs['cartoon_colors'],
                bottom_space=-0.78,
                show_labels=False
            )

    # title
    ax.set_title(
        temp_kwargs['title'],
        fontsize=12,
        fontname='Arial',
        color='k',
        pad=title_pad
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()
