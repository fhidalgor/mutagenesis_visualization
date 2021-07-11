#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


import numpy as np
import seaborn as sns
import pandas as pd
import copy
import matplotlib.pyplot as plt
from os import path
import os
from pathlib import Path
from typing import Union

# local modules
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils


# # Plot Functions

# ## Kernel

# In[41]:


def plot_kernel(
    self,
    cumulative=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    Plot univariate or bivariate distributions using kernel density estimation.

    Parameters
    ----------
    self : object from class *Screen*

    cumulative : bool, optional, default False
        If True, estimate a cumulative distribution function.
    
    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).

    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object=True. By default they do
        not get returned.

    """
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
    temp_kwargs['x_label'] = kwargs.get('x_label', r'$∆E^i_x$')
    temp_kwargs['y_label'] = kwargs.get('y_label', 'Probability density')

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # import parameters
    code_kwargs._parameters()

    # plot
    ax = sns.kdeplot(
        self.dataframe['Score_NaN'],
        cumulative=cumulative,
        color='red',
        lw=2,
    )

    # tune graph
    plt.xlabel(
        temp_kwargs['x_label'],
        fontsize=10,
        fontname='Arial',
        color='k',
        labelpad=0
    )
    plt.ylabel(['y_label'],
               fontsize=10,
               fontname='Arial',
               color='k',
               labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()


# ## Multiple Kernel plots

# In[47]:


def plot_multiplekernel(
    dict_entries,
    colors=['k', 'crimson', 'dodgerblue', 'g', 'silver'],
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a kernel density plot for multiple objects passed as a dictionary.
    If specified it can also draw a histogram. Uses sns.dispplot. Can manage either 
    Screen objects or dataframes out of the calculate_enrichments function.

    Parameters
    ----------
    dict_entries : dictionary containing dataframes
        Allows for either putting multiple objects as inputs or to use dataframes 
        that come out of the calculate_enrichments function. If you use an object, 
        you need to say object.dataframe.

    colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
        List of the colors (in order of arguments) that the kernels will have.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
        
    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
            
    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object=True. By default they do
        not get returned.    
        
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))

    # hard copy of data
    dict_copy = copy.deepcopy(dict_entries)

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # import parameters
    code_kwargs._parameters()

    # plot (allows two types of data input)
    for (label, dataset, color) in zip(dict_copy.keys(), dict_copy.values(),
                                       colors[0:len(dict_copy)]):
        if isinstance(dataset,
                      pd.core.frame.DataFrame):  # check if input is a dataframe
            # plot objects scores
            ax = sns.kdeplot(
                dataset['Score_NaN'], color=color, lw=2, label=label
            )
        else:  # otherwise assume its an array
            # get rid of stop codons
            dataset.drop('*', errors='ignore', inplace=True)
            dataset = dataset.stack()
            # plot stacked matrix
            ax = sns.kdeplot(
                dataset[~np.isnan(dataset)], color=color, lw=2, label=label
            )

    # tune graph
    plt.xlabel(
        r'$∆E^i_x$', fontsize=10, fontname='Arial', color='k', labelpad=0
    )
    plt.ylabel(
        'Probability density',
        fontsize=10,
        fontname='Arial',
        color='k',
        labelpad=3
    )
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()
    plt.legend(
        dict_copy.keys(),
        loc='best',
        frameon=False,
        fontsize=9,
        handlelength=1,
        handletextpad=0.5
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()


# ## Plot Histogram

# In[ ]:


def plot_hist(
    self,
    population='All',
    loc='upper left',
    output_file: Union[None, str, Path] = None,
    **kwargs
):
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
        r'$∆E^i_x$'
        if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.ylabel(
        'Probability density',
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=3
    )
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k')

    # axes limits. spacer will be 1 or the
    plt.xlim(temp_kwargs['xscale'])
    plt.xticks(
        np.arange(
            temp_kwargs['xscale'][0],
            temp_kwargs['xscale'][1] + temp_kwargs['tick_spacing'],
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

