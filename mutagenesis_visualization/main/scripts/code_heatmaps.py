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
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from os import path
from pathlib import Path
from typing import Union
from collections import defaultdict, Counter
from scipy.cluster import hierarchy

# local modules
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils


# # Plot Functions

# ## Heatmap

# ### Full Heatmap

# In[7]:


def plot_heatmap(
    self,
    nancolor='lime',
    show_cartoon=False,
    show_snv=False,
    hierarchical=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a heatmap plot of the enrichment scores.

    Parameters
    ----------
    self : object from class *Screen*
    
    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.
    
    show_carton : boolean, default False
        If true, the plot will display a cartoon with the secondary structure. 
        The user must have added the secondary structure to the object. 
    
    show_snv : boolean, default False
        If true, it will only display mutants that are a single nucleotide variant (SNV) of the wild-type
        protein sequence. The algorithm does not take into account the wild-type DNA allele, so it 
        will include any possible mutant that is one base away.
    
    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
                
    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).        
    Returns
    ----------
    fig, ax, ax2, ax3, averageresidue : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    '''
    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # sort data by rows in specified order by user
    df = code_utils._df_rearrange(
        code_utils._add_SNV_boolean(self.dataframe_stopcodons.copy()),
        temp_kwargs['neworder_aminoacids'],
        values='Score_NaN',
        show_snv=show_snv
    )

    # average of residues by positon
    average = [
        code_utils._add_SNV_boolean(self.dataframe.copy()
                                    ).groupby(by='Position').mean()['Score_NaN']
    ]

    # Create new sequence that we may to change order later
    self.sequence_updated = self.sequence

    # sort data by columns if specified by user
    if hierarchical:
        sorted_columns = _hierarchical_sort(df)
        sorted_columns_corrected = sorted_columns + self.start_position  # adds the initial number so index match
        df = df[sorted_columns_corrected]
        temp = pd.DataFrame(
            average
        )  # had to convert back to df to be able to sort properly
        average = temp[sorted_columns_corrected]
        # Change the sequence order
        self.sequence_updated = pd.DataFrame(list(self.sequence)
                                             ).T[sorted_columns]
        self.sequence_updated = list(self.sequence_updated.iloc[0])

    # declare figure and subplots
    figwidth = 14 * len(df.columns) / 165

    # Change parameters depending on whether cartoon is on or off
    if show_cartoon:
        figheight = 2.45
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = gridspec.GridSpec(
            nrows=3,
            ncols=2,
            height_ratios=[len(df), 1, 5],
            width_ratios=[len(df.columns), 1]
        )
    else:
        figheight = 2
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = gridspec.GridSpec(
            nrows=2,
            ncols=2,
            height_ratios=[len(df), 1],
            width_ratios=[len(df.columns), 1]
        )

    ax = plt.subplot(gs[0, 0])
    averageresidue = plt.subplot(gs[1, 0])
    cbar1 = plt.subplot(gs[0, 1])
    cbar2 = plt.subplot(gs[1, 1])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax.pcolormesh(
        df,
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=cmap,
        edgecolors='k',
        linewidths=0.2,
        antialiased=True,
        color='darkgrey'
    )

    # average by position
    heatmapaverageresidues = averageresidue.pcolormesh(
        average,
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=cmap,
        edgecolors='k',
        linewidths=0.2,
        antialiased=True,
        color='darkgrey'
    )

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False)
    ax.set_yticks(np.arange(len(df)) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # make new axes
    ax2 = ax.twiny()
    ax3 = ax.twinx()

    # tune the axes
    ax2.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False)
    ax3.set_yticks(np.arange(len(df)) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)
    ax3.tick_params(direction='out', pad=0.4)
    averageresidue.tick_params(direction='out', pad=-2)
    averageresidue.set_xticks(
        np.arange(len(df.columns)) + 0.5,
        minor=False,
    )
    averageresidue.set_yticks(np.arange(0.5) + 0.5)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())
    ax3.set_ylim(ax.get_ylim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax3.invert_yaxis()

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')
    ax3.yaxis.set_ticks_position('none')
    averageresidue.xaxis.set_ticks_position('none')
    averageresidue.yaxis.set_ticks_position('none')

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(
        list(self.sequence_updated),
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False
    )
    ax.set_yticklabels(
        temp_kwargs['neworder_aminoacids'],
        fontsize=6,
        fontname="Arial",
        color='k',
        minor=False
    )
    # For numbering labels, change if hierarchical sorting is true
    if not (hierarchical):
        ax2.set_xticklabels(
            temp_kwargs['number_sequencelabels'][0:len(df.columns)],
            fontsize=10,
            fontname="Arial",
            color='k',
            minor=False
        )
    else:
        ax2.tick_params(direction='out', pad=7)
        ax2.set_xticklabels(
            sorted_columns_corrected,
            fontsize=5,
            fontname="Arial",
            color='k',
            minor=False,
            rotation=90
        )
    ax3.set_yticklabels(
        temp_kwargs['neworder_aminoacids'],
        fontsize=6,
        fontname="Arial",
        color='k',
        minor=False
    )
    averageresidue.set_xticklabels(
        list(self.sequence_updated),
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False
    )
    rowaverage = ''
    averageresidue.set_yticklabels(
        rowaverage, fontsize=6, fontname="Arial", color='k', minor=False
    )

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')
    for ylabel in ax3.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # for coloring the residues that are 10,20...
    for xtick, color in zip(ax.get_xticklabels(),
                            temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)
    for xtick, color in zip(averageresidue.get_xticklabels(),
                            temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)
    # _____________________________________________________________________________

    # for color bar format
    cbar1.axis('off')
    cbar2.axis('off')
    cb = plt.colorbar(
        heatmap,
        fraction=1,
        pad=0,
        ax=[cbar1],
        aspect=5,
        ticks=[
            temp_kwargs['colorbar_scale'][0],
            np.mean(temp_kwargs['colorbar_scale']),
            temp_kwargs['colorbar_scale'][1]
        ],
        orientation='vertical'
    )
    cb.ax.set_yticklabels(
        cb.ax.get_yticklabels(), fontsize=6, fontname="Arial", color='k'
    )
    cb.update_ticks()
    plt.text(
        1.2 + 10 / len(df.columns),
        0.7,
        r'$\langle∆E^x_i\rangle_x$',
        transform=cbar1.transAxes,
        horizontalalignment='center',
        fontsize=7,
        fontname="Arial",
        color='k'
    )

    gs.update(hspace=0.1, wspace=0.1 / len(df.columns) * 50)

    # for putting title on graph
    plt.title(
        temp_kwargs['title'],
        horizontalalignment='center',
        fontname="Arial",
        fontsize=12
    )

    # Cartoon
    if show_cartoon:
        _generate_cartoon(self, gs, 2, temp_kwargs['cartoon_colors'], 0.025)

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, ax2, ax3, averageresidue

    if temp_kwargs['show']:
        plt.show()


# In[3]:


def _labels(start_position=1):
    # residue label and color
    emptylist = [''] * 1000
    number_sequencelabels = list([
        'b' if index in np.arange(10 - (start_position % 10), 1000, 10) else 'k'
        for index, x in enumerate(emptylist)
    ])
    color_sequencelabels = list([
        index + start_position
        if index in np.arange(10 - (start_position % 10), 1000, 10) else ''
        for index, x in enumerate(emptylist)
    ])
    return number_sequencelabels, color_sequencelabels


def _generate_cartoon(
    self,
    gs,
    n_row,
    colors,
    bottom_space=0,
    fig_inches=13.91,
    show_labels=True
):
    '''Generates cartoon for heatmap'''
    # Create subplot
    cartoon = plt.subplot(gs[n_row, 0])

    # Generate coordinates of labels
    labels = list(Counter(self.secondary).keys())
    length = list(Counter(self.secondary).values())
    cumsum = length[:-1]
    cumsum.insert(0, self.start_position)
    cumsum = np.cumsum(cumsum)

    # Create cartoon
    for label, length, cum in zip(labels, length, cumsum):
        if 'β' in label:
            loopstructure = _loop(cum, length, color=colors[2])
            cartoon.add_patch(loopstructure)
            sheetstructure = _sheet(cum, length, colors[0])
            cartoon.add_patch(sheetstructure)
            x_label = cum + length - 3.5
            if length > 2 and show_labels:  # If beta sheet is too small, label does not fit
                if length == 3:
                    cartoon.text((x_label + 0.6),
                                 -0.25,
                                 label,
                                 name='Arial',
                                 fontweight='normal',
                                 size=8.5 * fig_inches / 13.91,
                                 multialignment='right')
                else:
                    cartoon.text((x_label),
                                 -0.25,
                                 label,
                                 name='Arial',
                                 fontweight='normal',
                                 size=8.5 * fig_inches / 13.91,
                                 multialignment='right')
        elif 'α' in label:
            helixstructure = _helix(cum, length, colors[1])
            cartoon.add_patch(helixstructure)
            x_label = cum + length / 2 - 1
            if length > 2 and show_labels:
                cartoon.text((x_label),
                             -0.3,
                             label,
                             name='Arial',
                             fontweight='normal',
                             size=9 * fig_inches / 14,
                             multialignment='center')
        elif 'L' in label:
            loopstructure = _loop(cum, length, colors[2])
            cartoon.add_patch(loopstructure)

    # format of secondary cartoon
    cartoon.xaxis.set_ticks_position('none')
    cartoon.yaxis.set_ticks_position('none')
    cartoon.axis('off')

    # size
    cartoon.set_xlim(
        self.start_position - 0.1,
        len(self.secondary) + self.start_position + 0.2
    )
    cartoon.set_ylim(-2, 2.5)

    # adjust proximity to heatmap
    box = cartoon.get_position()
    box.y0 = box.y0 - bottom_space
    box.y1 = box.y1 - bottom_space
    cartoon.set_position(box)

    return


def _sheet(starting_aa, length_aa, color='lightgreen'):
    dx = length_aa
    sheetstructure = patches.FancyArrow(
        starting_aa,
        0.25,
        dx,
        0,
        width=2,
        length_includes_head=True,
        head_width=4,
        head_length=3,
        shape='full',
        overhang=0,
        head_starts_at_zero=False,
        ec='k',
        fc=color
    )
    return sheetstructure


def _helix(starting_aa, length_aa, color='lavender'):
    """produces matplotlib.Rectangle of length specified by length_aa """
    dx = length_aa  # so i can overlap tip
    helixstructure = plt.Rectangle((starting_aa, -0.85),
                                   dx,
                                   2.2,
                                   fc=color,
                                   ec='k')
    return helixstructure


def _loop(starting_aa, length_aa, color='k'):
    dx = length_aa
    loopstructure = plt.Rectangle((starting_aa, 0), dx, 0.5, fc=color)
    return loopstructure


def _hierarchical_sort(df):
    """
    sorts columns of dataset using hierarchical clustering
    returns order to rearrange dataset
    
    """

    # replaces NaN values with 0
    df_1 = df.fillna(0)

    # give sorted column order
    Z = hierarchy.ward(df_1.T)
    new_order = hierarchy.leaves_list(Z)

    return new_order


# ### Grouped Heatmap

# In[4]:


def plot_heatmap_rows(
    self,
    selection=['E', 'Q', 'A', 'P', 'V', 'Y'],
    nancolor='lime',
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a heatmap plot enrichment scores of selected aminoacids. 

    Parameters
    ----------
    self : object from class *Screen*
    
    selection : list of aa to show, default ['E','Q','A','P','V','Y']. 
    
    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.
    
    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
                    
    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
    Returns
    ----------
    fig, ax, ax2 : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    '''
    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # Add group and pivot df
    df = code_utils._select_aa(
        self.dataframe_stopcodons, selection, values='Score_NaN'
    )
    dataset = df.to_numpy()

    # The size can be changed. I found it empirically
    figwidth = 14 * len(dataset[0]) / 165
    figheight = 2 / 21 * len(selection)
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[len(dataset[0]), 1])
    ax = plt.subplot(gs[0, 0])
    cbar1 = plt.subplot(gs[0, 1])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax.pcolormesh(
        dataset,
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=cmap,
        edgecolors='k',
        linewidths=0.2,
        antialiased=True,
        color='darkgrey'
    )

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # second axis
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(
        list(self.sequence),
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False
    )
    ax.set_yticklabels(
        list(df.T.columns),
        fontsize=6,
        fontname="Arial",
        color='k',
        minor=False
    )
    ax2.set_xticklabels(
        temp_kwargs['number_sequencelabels'][0:len(dataset[0])],
        fontsize=10,
        fontname="Arial",
        color='k',
        minor=False
    )

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # for coloring the residues that are 10,20...
    for xtick, color in zip(ax.get_xticklabels(),
                            temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)

    # for putting title on graph
    plt.title(
        temp_kwargs['title'],
        horizontalalignment='center',
        fontname="Arial",
        fontsize=12
    )

    # for color bar format
    cbar1.axis('off')
    cb = plt.colorbar(
        heatmap,
        fraction=1,
        pad=0,
        ax=[cbar1],
        aspect=5,
        ticks=[
            temp_kwargs['colorbar_scale'][0],
            np.mean(temp_kwargs['colorbar_scale']),
            temp_kwargs['colorbar_scale'][1]
        ],
        orientation='vertical'
    )
    cb.ax.set_yticklabels(
        cb.ax.get_yticklabels(), fontsize=8, fontname="Arial", color='k'
    )
    cb.update_ticks()
    gs.update(hspace=0.1, wspace=0.1 / len(dataset[0]) * 50)

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, ax2

    if temp_kwargs['show']:
        plt.show()


# ### Subset Heatmap

# In[5]:


def plot_heatmap_columns(
    self,
    segment,
    ylabel_color='k',
    nancolor='lime',
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a heatmap plot enrichment scores but only plots a selected segment.

    Parameters
    ----------
    self : object from class *Screen*
    
    segment : list
        Segment is typed as [20,40] and includes both residues 20 and 40.
    
    ylabel_color : str, default 'k'
        Choose white if you don't want amino acid y axis label.
    
    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.
    
    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
                    
    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
        
    Returns
    ----------
    fig, ax, ax2 : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    '''

    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # sort data in specified order by user
    df_whole = code_utils._df_rearrange(
        self.dataframe_stopcodons,
        temp_kwargs['neworder_aminoacids'],
        values='Score_NaN'
    )

    # select subset
    c0 = segment[0] - self.start_position
    c1 = segment[1] - self.start_position + 1
    df = df_whole.iloc[:, c0:c1]

    # the size can be changed
    figwidth = 2 * len(df.columns) / 22
    figheight = 2
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    # needed to set autoscale off to avoid missalignment
    ax = plt.subplot(gs[0])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax.pcolormesh(
        df,
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=cmap,
        edgecolors='k',
        linewidths=0.2,
        antialiased=True,
        color='darkgrey'
    )

    # put the major ticks at the middle of each cell
    ax.set_xticks(
        np.arange(len(df.columns)) + 0.5,
        minor=False,
    )
    ax.set_yticks(np.arange(len(df)) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # second axis
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(
        list(self.sequence)[segment[0] - self.start_position:segment[1] -
                            self.start_position + 1],
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False
    )
    ax.set_yticklabels(
        temp_kwargs['neworder_aminoacids'],
        fontsize=6,
        fontname="Arial",
        color=ylabel_color,
        minor=False
    )

    ax2_label = (segment[1] - segment[0] + 1) * ['']
    ax2_label[0] = segment[0]
    ax2_label[-1] = segment[1]
    ax2.set_xticklabels(
        ax2_label, fontsize=7, fontname="Arial", color='k', minor=False
    )

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, ax2

    if temp_kwargs['show']:
        plt.show()

