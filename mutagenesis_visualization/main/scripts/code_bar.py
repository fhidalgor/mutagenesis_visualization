#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


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

# Local imports
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
    from mutagenesis_visualization.main.scripts.code_heatmaps import _generate_cartoon
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils
    from code_heatmaps import _generate_cartoon


# # Plot Functions

# ## Mean Plots

# ### Bar graph Enrichment

# In[7]:


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

    show_carton : boolean, default False
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


def _parameters_mean():
    # normal font
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'

    # Parameters for all graphs
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9

    return


# ### Compare two proteins

# In[6]:


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

    show_carton : boolean, default False
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


# ### Bar graph Counts

# In[9]:


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


# ### Positional

# In[10]:


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


# ## Stacked bar cumulative %AA representation

# In[11]:


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


# In[12]:


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

