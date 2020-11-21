#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[ ]:


import numpy as np
import seaborn as sns
import pandas as pd
import copy
import matplotlib.pyplot as plt
from os import path
import os
from pathlib import Path
from typing import Union

# kwargs and parameters
try:
    import import_notebook
except ModuleNotFoundError:
    import sys
    sys.path.append('mutagenesis_visualization/main/scripts/')

import code_kwargs
import code_utils


# # Plot Functions

# ## Kernel

# In[ ]:


def plot_kernel(self, kernel='gau', kernel_label='KDE', histogram=False,
                fit=None, fit_label='_nolegend_', extra_dist=None,
                extra_dist_label='_nolegend_',  
                output_file: Union[None, str, Path] = None,**kwargs):
    '''
    Generate a kernel density plot. If specified it can also draw a histogram. Uses sns.distplot.

    Parameters
    ----------
    self : object from class *Screen*
    
    kernel : str, default gau
        options are ['biw','cos','epa','gau','tri','triw']
    
    kernel_label : str, default '_nolegend_'
    
    histogram : boolean, default False
    
    fit : boolean, optional
        ask sns.distplot to fit a function
    
    fit_label : str, default '_nolegend_'
    
    extra_dist : [x,y], optional
        fit any distribution you want. Input the x and y coordinates
    
    extra_dist_label : str, default '_nolegend_'
    
    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
        
    **kwargs : other keyword arguments
        
    Returns
    ----------
    None
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # import parameters
    code_kwargs._parameters()

    # plot
    ax = sns.distplot(self.dataframe['Score_NaN'], kde=True, hist=histogram, norm_hist=True,
                      kde_kws={'kernel': kernel, 'color': temp_kwargs['color'], 'lw': 2, 'label': kernel_label}, fit=fit,
                      fit_kws={'label': fit_label, 'linestyle': 'dotted', 'color': 'red'})

    # plot extra distribution
    if extra_dist is not None:
        plt.plot(extra_dist[0], extra_dist[1], linewidth=2, linestyle='dotted',
                 color='green', label=extra_dist_label)

    # tune graph
    plt.xlabel(r'$∆E^i_x$', fontsize=10,
               fontname='Arial', color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,
               fontname='Arial', color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()
    ax.legend(loc='best', frameon=False, fontsize=9,
              handlelength=1, handletextpad=0.5)

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()

    return


# ## Multiple Kernel plots

# In[ ]:


def plot_multiplekernel(dict_entries, kernel='gau',
                        colors=['k', 'crimson', 'dodgerblue', 'g', 'silver'], 
                        output_file: Union[None, str, Path] = None,**kwargs):
    '''
    Generate a kernel density plot for multiple objects passed as a dictionary.
    If specified it can also draw a histogram. Uses sns.distplot. Can manage either Screen objects
    or dataframes out of the calculate_enrichments function.

    Parameters
    ----------
    dict_entries : dictionary containing dataframes
        Allows for either putting multiple objects as inputs or to use dataframes 
        that come out of the calculate_enrichments function. If you use an object, 
        you need to say object.dataframe.

    kernel : str, default gau
        options are ['biw','cos','epa','gau','tri','triw'].

    colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
        List of the colors (in order of arguments) that the kernels will have.

output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
        
    **kwargs : other keyword arguments

    Returns
    ----------
    None
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
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
    for (label, dataset, color) in zip(dict_copy.keys(), dict_copy.values(), colors[0:len(dict_copy)]):
        if 'Score' in dataset.columns:
            # plot objects scores
            sns.distplot(dataset['Score_NaN'], hist=False,
                         kde_kws={"color": color, "lw": 2, "label": label})
        else:
            # get rid of stop codons
            dataset.drop('*', errors='ignore', inplace=True)
            dataset = dataset.stack()
            # plot stacked matrix
            sns.distplot(dataset[~np.isnan(dataset)], kde=True, hist=False, 
                         kde_kws={"color": color, "lw": 2, "label": label})

    # tune graph
    plt.xlabel(r'$∆E^i_x$', fontsize=10,
               fontname='Arial', color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,
               fontname='Arial', color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()
    plt.legend(dict_copy.keys(), loc='best', frameon=False, fontsize=9,
               handlelength=1, handletextpad=0.5)

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()

    return


# ## Plot Histogram

# In[ ]:


def plot_hist(self, population='All', loc='upper left', 
              output_file: Union[None, str, Path] = None,**kwargs):
    '''
    Generate a histogram plot. Can plot single nucleotide variants (SNVs) or non-SNVs only

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
        bins : int, default 50. 
            Number of bins for the histogram.
    Returns
    ----------
    None.

    '''
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
    plt.xlabel(r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
               fontsize=10, fontname="Arial", color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,
               fontname="Arial", color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k')

    # axes limits. spacer will be 1 or the
    plt.xlim(temp_kwargs['xscale'])
    plt.xticks(np.arange(temp_kwargs['xscale'][0], temp_kwargs['xscale']
                         [1]+temp_kwargs['tick_spacing'], temp_kwargs['tick_spacing']))
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)
    if temp_kwargs['show']:
        plt.show()
    return

