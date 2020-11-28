#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


import numpy as np
import seaborn as sns
import pandas as pd
import itertools
import copy
from scipy import stats
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from sklearn import metrics
import freesasa
from os import path
import os
from pathlib import Path
from typing import Union

# local modules
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
    from mutagenesis_visualization.main.scripts.code_heatmaps import _labels
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils
    from code_heatmaps import _labels


# # Plot Functions

# ## Rank 

# In[2]:


def plot_rank(
    self,
    mode='pointmutant',
    outdf=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    DEPRECATED.
    Generate a rank plot so every mutation/residue is sorted based
    on enrichment score.

    Parameters
    ----------
    self : object from class *Screen*

    mode : str, default 'pointmutant'.
        Alternative set to "mean" for the mean of each position

    outdf : boolean, default False
        If set to true, will return the df with the rank of mutations

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
    
    df : Pandas dataframe
        Contains the mutations sorted by their performance.
        Needs to have outdf==True. By default they do
        not get returned.
    
    """
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (4, 2))
    temp_kwargs['x_label'] = kwargs.get('x_label', 'Rank')
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # Sort by enrichment scores
    df = self.dataframe.sort_values(by=['Score']).copy()

    # Chose mode:
    if mode == 'mean':
        df = df.groupby(by=['Position'], as_index=False).mean()
        df.sort_values(by=['Score'], inplace=True)

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # import parameters
    code_kwargs._parameters()

    # Scatter data points
    plt.scatter(np.arange(len(df), 0, -1), df['Score'], c='k', s=1)

    # Titles
    plt.title(
        temp_kwargs['title'], fontsize=12, fontname='Arial', color='k', pad=8
    )
    # Labels
    plt.ylabel(
        temp_kwargs['y_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k')

    # other graph parameters
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()

    if outdf:
        return df


# ## Miniheatmap

# ### Mean substitution heatmap

# In[3]:


def plot_miniheatmap(
    self,
    offset=0,
    background_correction=True,
    output_file: Union[None, str, Path] = None,
    **kwargs,
):
    """
    Generate a miniheatmap plot enrichment scores of mutagenesis selection
    assays.

    Parameters
    ----------
    self : object from class *Screen*

    offset : int, default 0
        Will group columns by residues. If the offset is not 0, it will use the values
        of the n+offset to group by. For example, you may want to see what happens when
        you have a Proline in front of the mutated residue. The algorithm can report
        the difference between the calculated value and the mean score for that particular
        substitution.
        Offset of 1 means that you evaluate the effect of following residue n+1 on n.
        Offset of -1 means that you look at the previous residue (n-1 on n).

    background_correction : boolean, default True
        If offset is nonzero, whether subtract the average effect of a substitution or not.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
            
    Returns
    ----------
    fig, ax, cb : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    """

    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # do offset if appropriate
    dataframe_stopcodons = _transform_dataset_offset(self, offset)

    # calculate condensed heatmap
    dataset = _condense_heatmap(
        dataframe_stopcodons, temp_kwargs['neworder_aminoacids']
    )

    # if offset is not 0
    if background_correction and offset != 0:
        if '*' in temp_kwargs['neworder_aminoacids']:
            temp_kwargs['neworder_aminoacids'].remove('*')
        # do offset, no stop codons
        dataset = _normalize_neighboreffect(
            self, offset, temp_kwargs['neworder_aminoacids']
        )

    fig, ax, cb = _plot_miniheatmap(dataset, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, cb


def _condense_heatmap(df, new_order):
    '''
    Converts the np.array with stored enrichment scores into the condensed heatmap
    '''
    # Convert dataset to df
    df = df.copy()
    df.drop(['Position'], axis=1, inplace=True)

    # Group by sequence and aminoacid, and then pivot table
    df_grouped = df.groupby(['Sequence', 'Aminoacid'], sort=False).mean()
    df_pivoted = df_grouped.pivot_table(
        values='Score', index='Aminoacid', columns='Sequence'
    )
    df_pivoted.reset_index(drop=False, inplace=True)

    # Sort in y axis desired order
    df_pivoted['Aminoacid'] = pd.Categorical(df_pivoted['Aminoacid'], new_order)
    df_pivoted = df_pivoted.sort_values(by=['Aminoacid'])

    # Sort in x axis desired order
    x_order = code_utils._common(new_order, list(df_pivoted.columns))

    # Drop amino acid column
    data_dropped = df_pivoted.drop(['Aminoacid'], axis=1)

    return data_dropped[x_order]


def _offset_sequence(dataset, sequence, start_position, offset):
    """
    Internal function that offsets the input sequence.

    Parameters
    -----------
    dataset, sequence, start_position, offset

    Returns
    --------
    string containing trimmed sequence

    """
    # Deep copy sequence
    sequence = copy.deepcopy(sequence)

    # truncate sequence
    if offset > 0:
        sequence = sequence + 'X' * np.absolute(offset)
        trimmedsequence = sequence[start_position - 1 + offset:len(dataset[0]) +
                                   start_position - 1 + offset]
    else:
        sequence = 'X' * (np.absolute(offset)) + sequence
        trimmedsequence = sequence[start_position - 1:len(dataset[0]) +
                                   start_position - 1]

    return trimmedsequence


def _transform_dataset_offset(self, offset, stopcodons=True):
    '''
    Generate a dataframe with the sequence offset. Reutilizes _transform_dataset
    '''
    # Add offset sequence
    offset_sequence = _offset_sequence(
        self.dataset, self.sequence_raw, self.start_position, offset
    )
    df = self.dataframe_stopcodons.copy(
    ) if stopcodons is True else self.dataframe.copy()

    # Copy old sequence
    df['Sequence_old'] = df['Sequence']
    # Count amino acids
    aa_number = len(set(df['Aminoacid']))
    # Generate new offset sequence
    df['Sequence'] = np.ravel([[aa] * aa_number for aa in offset_sequence])

    # Drop rows with X
    df.drop(df.index[df['Sequence'] == 'X'], inplace=True)

    return df


# ### Neighbor residues

# In[4]:


def plot_neighboreffect (self, offset=1, output_file: Union[None, str, Path] = None,
                         **kwargs):
   """
   DEPRECATED.

   Generate a miniheatmap plot telling you the effect of having a residue in front or behind.
   It corrects for the effect of that amino acid on the rest of the population.

   Parameters
   ----------
   self : object from class *Screen*

   offset : int, default 1
      if you want to study effects of a residue when is behind or in front of another residue.
      offset of 1 means that you evaluate the effect of following residue n+1 on n. On a "MTEY..." sequence,
      you would look at the effect of T on M, E on T, Y on E, etc.. and then group by residue (n+1).
      offset of -1 means that you look at the previous residue (n-1 on n).

   output_file : str, default None
      If you want to export the generated graph, add the path and name of the file.
      Example: 'path/filename.png' or 'path/filename.svg'.

   **kwargs : other keyword arguments
       return_plot_object : boolean, default False
           If true, will return plotting objects (ie. fig, ax).
           
   Returns
   ----------
   fig, ax, cb : matplotlib figure and subplots
       Needs to have return_plot_object==True. By default they do
       not get returned.
       
   """
   # load font parameters
   code_kwargs._font_parameters()

   # update kwargs
   temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
   temp_kwargs.update(kwargs)
   if '*' in temp_kwargs['neworder_aminoacids']: temp_kwargs['neworder_aminoacids'].remove('*')
   
   # do offset, no stop codons
   df = _normalize_neighboreffect(self,offset,temp_kwargs['neworder_aminoacids'])
   
   # Plot
   fig, ax, cb = _plot_miniheatmap(df,temp_kwargs)
   
   # return matplotlib object
   if temp_kwargs['return_plot_object']:
       return fig, ax, cb
   
   # show figure    
   if temp_kwargs['show']: 
       plt.show()
   
def _plot_miniheatmap(df,output_file, temp_kwargs):
   """
   Aux plot that will do the heavy lifting to plot a miniheatmap.
   
   Parameters
   ------------
   df: pandas dataframe
       dataframe with the data to plot.
   
   output_file : str, default None
       If you want to export the generated graph, add the path and name of the file.
       Example: 'path/filename.png' or 'path/filename.svg'.
   
   temp_kwargs : kwargs
   
   Returns
   --------
   fig, ax, cb : matplotlib figure and subplots
       
   """
   # declare figure and subplots
   coeff = len(df.columns)/19*1.05
   fig = plt.figure(figsize=(2.5*coeff, 2.5))
   gs = gridspec.GridSpec(nrows=1, ncols=1)
   ax = plt.subplot(gs[0])

   # main heatmap
   heatmap = ax.pcolor(df.to_numpy(), vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                       cmap=temp_kwargs['colormap'], edgecolors='k', linewidths=0.2, color='darkgrey')

   # ____________axes manipulation____________________________________________
   # put the major ticks at the middle of each cell
   ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
   ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

   # position of axis labels
   ax.tick_params('x', direction='out', pad=-2.5)
   ax.tick_params('y', direction='out', pad=0.4)

   # want a more natural, table-like display
   ax.invert_yaxis()
   ax.xaxis.tick_top()

   # remove ticks
   ax.xaxis.set_ticks_position('none')
   ax.yaxis.set_ticks_position('none')

   # so labels of x and y do not show up and my labels show up instead
   ax.set_xticklabels(list(df.columns), fontsize=6.5,
                      fontname="Arial", color='k', minor=False)
   ax.set_yticklabels(temp_kwargs['neworder_aminoacids'],
                      fontsize=6.5, fontname="Arial", color='k', minor=False)

   # align the labels of the y axis
   for ylabel in ax.get_yticklabels():
       ylabel.set_horizontalalignment('center')

   # _____________________________________________________________________________

   # for color bar format
   cb = plt.colorbar(heatmap, fraction=0.025, pad=0.05, aspect=5, ticks=[temp_kwargs['colorbar_scale'][0], np.mean(
       temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]], orientation='vertical')
   cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k')
   cb.update_ticks()
   plt.text(len(df.columns)+2, 7.8, r'$\langle∆E^x_i\rangle_x$', horizontalalignment='center',
            fontsize=7, fontname="Arial", color='k')

   # for putting title on graph
   plt.title(temp_kwargs['title'], horizontalalignment='center',
             fontname="Arial", fontsize=10, pad=10)
   plt.ylabel('Amino Acid Substitution', fontsize=10, labelpad=-1)

   # save file
   code_utils._save_work(fig, output_file, temp_kwargs)
   
   # return matplotlib object
   return fig, ax, cb


def _normalize_neighboreffect(self,offset,neworder):
   '''
   For every residue, subtract the average effect of a substitution
   Returns a normalized dataframe
   '''
   aalist = list('ACDEFGHIKLMNPQRSTVWY')
   # Add offset sequence to df
   df = _transform_dataset_offset(self,offset,False)
   
   # calculate mean effect using condensed heatmap
   mean = _condense_heatmap(self.dataframe, aalist)
   
   df_normalized = pd.DataFrame()
   for aa in aalist:
       # Choose the neighbors of an aa
       aa_neighbors = df.loc[df['Sequence']==aa]
       # Do the mean substitution of amino acids that are repeated
       aa_neighbors = aa_neighbors.groupby(['Sequence_old','Aminoacid'],as_index=False).mean()
       # Make into table
       aa_neighbors_pivoted = aa_neighbors.pivot_table(values='Score', index='Aminoacid',  columns='Sequence_old')
       aa_neighbors_pivoted.reset_index(drop=True, inplace=True)
       # Get the mean of the amino acids that appear in the aa_neighbors subset
       mean_neighbors = mean[list(aa_neighbors_pivoted.columns)]
       # Subtract average effect and do mean
       df_normalized[aa] = (aa_neighbors_pivoted - mean_neighbors).mean(axis=1)
   
   # Sort by aa
   df_normalized = df_normalized[neworder]
   # Sort in y axis desired order
   df_normalized = _sort_yaxis_aminoacids(df_normalized,neworder,aalist)
   return df_normalized

def _sort_yaxis_aminoacids(df,neworder,oldorder=list('ACDEFGHIKLMNPQRSTVWY')):
   # Sort in y axis desired order
   df['Aminoacid_new'] = oldorder
   df['Aminoacid_new'] = pd.Categorical(df['Aminoacid_new'], neworder)
   df.sort_values(by=['Aminoacid_new'],inplace=True)
   df.drop(['Aminoacid_new'], inplace=True, axis=1)
   
   return df


# ## Secondary Structure

# In[5]:


def plot_secondary(self, output_file: Union[None, str, Path] = None, **kwargs):
    """
    Generates a bar plot of data sorted by secondary elements (alpha helices
    and beta sheets).

    Parameters
    -----------
    self : object from class *Screen*

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
    # Load parameters
    code_kwargs._parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-2, 1))

    # Get data
    df = _calculate_secondary(self.dataframe, self.secondary_dup)

    # Color
    df['Color'] = df.apply(
        code_utils._color_data,
        axis=1,
        args=(temp_kwargs['color_gof'], temp_kwargs['color_lof'])
    )

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ticks = np.arange(0, len(df))  # label locations
    width = 0.5
    labels = df['Secondary']

    # Plot figure
    ax.bar(
        ticks,
        df['Score'],
        width,
        color=df['Color'],
        ec='k',
    )

    # graph parameters
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        labels,
        fontsize=9,
        fontname="Arial",
        color='k',
        minor=False,
        rotation=0
    )
    ax.set_ylabel(
        r'$∆E^i_x$',
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=12,
        rotation=0
    )
    ax.set_ylim(temp_kwargs['yscale'])
    plt.title(
        temp_kwargs['title'],
        horizontalalignment='center',
        fontname="Arial",
        fontsize=10,
        pad=5
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()


def _calculate_secondary(df, secondary):
    '''
    Returns copy
    '''
    df = df.copy()
    df.insert(4, 'Secondary', secondary)
    df = df.groupby(['Secondary'], as_index=False, sort=False).mean()
    df = df[df['Secondary'].str.startswith(('β', 'α'))]
    df = df.drop(['Position'], axis=1)
    return df


# ## ROC AUC

# In[6]:


def plot_roc(
    self, df_class=None, output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Generates ROC AUC plot. It compares enrichment scores to some labels that
    the user has specified.

    Parameters
    -----------
    self : object from class *Screen*

    df_class: Pandas dataframe
        A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
        containing the true class of that mutation. The true class can also be an input when creating the object.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    Returns
    --------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """
    # Use default class
    if df_class is None:
        df_class = self.roc_df

    # Merge dataframe with classes
    df = _mergeclassvariants(df_class, self.dataframe)

    # Calculate ROC parameters
    fpr, tpr, auc, _ = _rocauc(df)

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2.5))

    # import parameters
    code_kwargs._parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    lw = 2
    plt.plot(fpr, tpr, color='k', lw=lw, label='AUC = %0.2f' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

    # Graph limits
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    tick_spacing = 0.2
    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(tick_spacing)
    )  # Plt ticks
    ax.yaxis.set_major_locator(
        ticker.MultipleLocator(tick_spacing)
    )  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel(
        'True Positive Rate',
        fontsize=12,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.xlabel('False Positive Rate', fontsize=12, fontname="Arial", color='k')

    # Legend
    plt.legend(
        loc='lower right', handlelength=0, handletextpad=0, frameon=False
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()


def _rocauc(df):
    '''
    Calculate roc rates and auc.

    The input is a dataframe that contains [Variants,Class,Score]
    '''
    fpr, tpr, thresholds = metrics.roc_curve(
        df['Class'], df['Score'], drop_intermediate=True
    )
    auc = metrics.roc_auc_score(df['Class'], df['Score'])
    return fpr, tpr, auc, thresholds


def _mergeclassvariants(df_score, df_class):
    '''
    Merge the input dataframe containing the class (true score) for variants and the enrichment scores
    '''
    # Merge DMS with true score dataset
    df_merged = pd.merge(df_class, df_score, on=['Variant'], how='left')

    # Drop rows with Nan values
    df_merged.dropna(inplace=True)

    return df_merged


def _concattrueposneg(df_tp, df_tn, subset='Variant', keep='first'):
    '''
    Concat a df containing the true positive variants and the true negative variants

    Parameters
    -----------
    df_tp : Dataframe with the true positives
    df_tn : Dataframe with the true negatives
    subset : str, default Variant
    keep : {‘first’, ‘last’, False} 

    Returns
    --------
    None.
    '''
    # Concatenate tp and tn datasets
    df_true = pd.concat([df_tp, df_tn], sort=False)

    # Will keep a variant as true positive if found in both datasets (because could be a mistake in gnomAD)
    df_true.drop_duplicates(subset=subset, keep=keep, inplace=True)

    return df_true


# ## Cumulative

# In[7]:


def plot_cumulative(
    self, mode='all', output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Generates a cumulative plot of the enrichment scores by position.

    Parameters
    -----------
    self : object from class *Screen*

    mode : str, default 'all'
        Options are 'all','SNV' and 'nonSNV'.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
            
    Returns
    --------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    """

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2))
    temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 20)

    # import parameters
    code_kwargs._parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # Get data filtered
    df = _filter_bySNV(self, mode)
    cumsum = df.cumsum(skipna=False)['Score']
    plt.plot(df['Position'], cumsum / list(cumsum)[-1], color='red', lw=2)

    # y label
    y_label = 'Cumulative LoF'
    if list(cumsum)[-1] > 0:
        y_label = 'Cumulative GoF'

    # Graph limits
    plt.xlim(
        self.dataframe['Position'].min(), self.dataframe['Position'].max() + 1
    )
    plt.ylim(0, 1.1)

    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(temp_kwargs['tick_spacing'])
    )  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel(y_label, fontsize=12, fontname="Arial", color='k', labelpad=5)
    plt.xlabel('Position', fontsize=12, fontname="Arial", color='k', labelpad=0)

    # x=y line
    plt.plot([0, df['Position'].max()], [0, 1],
             color='silver',
             lw=2,
             linestyle='--')

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()


def _filter_bySNV(self, mode):

    # Select all, SNV, nonSNV
    if mode == 'all':
        df = self.dataframe
    elif mode == 'SNV':
        df = self.dataframe_SNV
    elif mode == 'nonSNV':
        df = self.dataframe_nonSNV
    df = df.groupby(by='Position', as_index=False).mean()
    return df


# ## Box Plot

# In[8]:


def plot_box(binned_x, y, output_file: Union[None, str, Path] = None, **kwargs):
    """
    Generates a boxplot. Data needs to be binned prior before using this
    function.

    Parameters
    -----------
    x, binned_y : arrays
        Contain the data is going to plot

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
    # Load parameters
    code_kwargs._parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # Plot data
    ax = sns.boxplot(binned_x, y, color='white', fliersize=2)

    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')

    # graph parameters
    plt.title(
        temp_kwargs['title'], fontsize=10, fontname='Arial', color='k', pad=8
    )
    plt.ylabel(
        temp_kwargs['y_label'],
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k')

    # axes limits
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()

