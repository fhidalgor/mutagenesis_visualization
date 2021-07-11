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


def plot_meancounts(self, output_file: Union[None, str, Path] = None, **kwargs):
    """
    Plot in a bargraph the mean counts for each residue of the protein.

    Parameters
    ----------
    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
        text_labels : list of lists, default empty
            If you want to add a label to the graph, add the coordinates and the text.
            Example: text_labels = [[x0,y0,text0],[x1,y1,text1]].

    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (10, 4))
    temp_kwargs['y_label'] = kwargs.get('y_label', 'Mean counts')

    # load parameters
    _parameters_mean()

    # make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ax = self.dataframe.mean().T.plot.bar(ax=ax, color='red')

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(
        temp_kwargs['y_label'],
        fontsize=12,
        fontname="Arial",
        color='k',
        labelpad=0,
        rotation=90
    )
    ax.set_xlabel(
        'Amino acid position',
        fontsize=12,
        fontname="Arial",
        color='k',
        labelpad=4
    )
    ax.set_xticklabels(self.positions)

    # Title
    plt.title(
        temp_kwargs['title'], fontsize=14, fontname='Arial', color='k', pad=0
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


def _inputtext(text_entries):
    '''the user can input text as a variable by manually giving the coordinates'''
    if text_entries:
        for entry in text_entries:
            plt.text(entry[0], entry[1], entry[2])
    return
