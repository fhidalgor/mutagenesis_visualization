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

def plot_library_representation(
    self, output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Generates a cumulative stacked bar plot. Each bar represents an amino acid
    position, and each color indicates the observed variant frequency.

    Parameters
    ----------
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
    temp_kwargs['figsize'] = kwargs.get('figsize', (10, 4))

    # load parameters
    _parameters_mean()

    # Transform data
    df_grouped = _group_codons_to_aa(self)
    df_percentage = _percentage_column(df_grouped)

    # colors
    colors = plt.cm.tab20(np.linspace(0, 1, 20))

    # plot
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ax = df_percentage.T.plot.bar(stacked=True, ax=ax, color=colors)

    # axes parameters
    ax.set_ylim(0, 100)
    ax.set_ylabel(
        'Cumulative % AA representation',
        fontsize=12,
        fontname="Arial",
        color='k',
        labelpad=0,
        rotation=90
    )
    ax.yaxis.set_major_formatter(ticker.PercentFormatter())

    ax.set_xlabel(
        'Amino acid position',
        fontsize=12,
        fontname="Arial",
        color='k',
        labelpad=4
    )
    ax.set_xticklabels(self.positions)

    plt.title(
        temp_kwargs['title'], fontsize=14, fontname='Arial', color='k', pad=20
    )

    # legend
    plt.legend(
        markerscale=0.5,
        frameon=False,
        framealpha=0,
        handlelength=0.75,
        handletextpad=0.25,
        ncol=len(df_percentage),
        columnspacing=0.75,
        bbox_to_anchor=(0.5, 1.09),
        loc='upper center'
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()


# In[ ]:


def _group_codons_to_aa(self):
    """
    Group different codons that are synonymous.

    Returns sum of counts

    """    # copy df
    df = self.dataframe.copy()

    df['Aminoacid'] = self.aminoacids

    # Group by mean
    df = df.groupby(as_index=True, by='Aminoacid', sort=False).sum()
    return df


def _percentage_column(df):
    """
    Make the percentage per column.
    """
    return df / df.sum(axis=0) * 100
