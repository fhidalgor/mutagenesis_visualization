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

def plot_mean(
    self,
    mode='mean',
    show_cartoon=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Plot in a bargraph the mean enrichment for each residue of the protein.
    Red for gain of function, blue for loss of function

    Parameters
    ----------
    self : object from class *Screen*

    mode : str, default 'mean'
        Specify what enrichment scores to show. If mode = 'mean', it will show the mean of
        each position. If mode = 'A', it will show the alanine substitution profile. Can be
        used for each amino acid. Use the one-letter code and upper case.

    show_cartoon : boolean, default False
        If true, the plot will display a cartoon with the secondary structure.
        The user must have added the secondary structure to the object.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).

        color_gof : str, default 'red'
            Choose color to color positions with an enrichment score > 0.

        color_lof : str, default 'blue'
            Choose color to color positions with an enrichment score < 0.

    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2.5))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # load parameters
    _parameters_mean()

    # Select grouping
    if mode == 'mean':
        df = self.dataframe.groupby('Position', as_index=False).mean()
    else:
        df = self.dataframe.loc[self.dataframe['Aminoacid'] == mode].copy()

    df['Color'] = df.apply(
        code_utils._color_data,
        axis=1,
        args=(temp_kwargs['color_gof'], temp_kwargs['color_lof'])
    )

    # make figure
    if show_cartoon:
        fig = plt.figure(figsize=temp_kwargs['figsize'])
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
        ax = plt.subplot(gs[0])
    else:
        fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    width = 1.2

    # Color based on values
    ax.bar(df['Position'], df['Score'], width, color=df['Color'], snap=False)

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
        _generate_cartoon(
            self,
            gs,
            1,
            temp_kwargs['cartoon_colors'],
            bottom_space=-0.78,
            show_labels=False
        )
        title_pad = 2.5

    # Plot title
    plt.title(
        temp_kwargs['title'],
        fontsize=12,
        fontname='Arial',
        color='k',
        pad=title_pad
    )

    # Put text labels
    _inputtext(temp_kwargs['text_labels'])

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()
