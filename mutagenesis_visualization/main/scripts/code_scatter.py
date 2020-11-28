#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[2]:


import numpy as np
import seaborn as sns
import pandas as pd
import itertools
import copy
from scipy import stats
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from os import path
import os
from pathlib import Path
from typing import Union
from collections import defaultdict

# local modules
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils


# # Plot Functions

# ## Scatter

# In[ ]:


def plot_scatter(
    self,
    obj2,
    mode='pointmutant',
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a scatter plot between object and a second object of the same class.

    Parameters
    ----------
    self : object from class *Screen*

    obj2 : object from class *Screen* to do the scatter with

    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position

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
        
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))

    # Chose mode:
    if mode == 'pointmutant':
        df = code_utils._process_bypointmutant(self, obj2)
    else:
        df = code_utils._process_meanresidue(self, obj2)

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # import parameters
    code_kwargs._parameters()

    # Scatter data points
    plt.scatter(
        df['dataset_1'],
        df['dataset_2'],
        c='k',
        s=8,
        alpha=0.5,
        rasterized=True,
        label='_nolegend_'
    )

    # Titles
    plt.title(
        temp_kwargs['title'], fontsize=12, fontname='Arial', color='k', pad=8
    )
    plt.ylabel(
        temp_kwargs['y_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k')

    # correlation and R2
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df['dataset_1'], df['dataset_2']
    )
    R2 = str(round(r_value**2, 2))
    legend_label = "$R^2$ = {}".format(R2)
    # fit and graph line
    fit = np.polyfit(df['dataset_1'], df['dataset_2'], 1)
    plt.plot(
        np.unique(df['dataset_1']),
        np.poly1d(fit)(np.unique(df['dataset_1'])),
        color='r',
        linewidth=1,
        label=legend_label
    )
    plt.grid()

    # other graph parameters
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(temp_kwargs['tick_spacing'])
    )
    ax.yaxis.set_major_locator(
        ticker.MultipleLocator(temp_kwargs['tick_spacing'])
    )
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()

    # Legend
    plt.legend(
        loc='upper left',
        handlelength=0,
        handletextpad=0,
        frameon=False,
        fontsize=10
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()

