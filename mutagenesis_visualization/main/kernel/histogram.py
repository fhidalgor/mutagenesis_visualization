import numpy as np
import seaborn as sns
import pandas as pd
import copy
import matplotlib.pyplot as plt
from os import path
import os
from pathlib import Path
from typing import Union


def plot_hist(self, population='All', loc='upper left', output_file: Union[None, str, Path] = None, **kwargs):
    """
    Generate a histogram plot. Can plot single nucleotide variants (SNVs) or
    non-SNVs only.

    Parameters
    ----------
    population : str, default 'All'.
        Other options are 'SNV' and 'nonSNV'.

    loc : str, default 'upper left'.
        Position of the legend.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
        bins : int, default 50.
            Number of bins for the histogram.

    Returns
    ----------
    fig : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))

    # Select case input data
    df = self.dataframe['Score_NaN']
    if population == 'SNV':
        df = self.dataframe_SNV['Score_NaN']
    elif population == 'nonSNV':
        df = self.dataframe_nonSNV['Score_NaN']

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # Import parameters
    code_kwargs._parameters()

    # plot figure
    plt.hist(df, density=True, bins=temp_kwargs['bins'], color='k')

    # axes labels and title
    plt.xlabel(
        r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.ylabel('Probability density', fontsize=10, fontname="Arial", color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k')

    # axes limits. spacer will be 1 or the
    plt.xlim(temp_kwargs['xscale'])
    plt.xticks(
        np.arange(
            temp_kwargs['xscale'][0], temp_kwargs['xscale'][1] + temp_kwargs['tick_spacing'],
            temp_kwargs['tick_spacing']
        )
    )
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig

    if temp_kwargs['show']:
        plt.show()
