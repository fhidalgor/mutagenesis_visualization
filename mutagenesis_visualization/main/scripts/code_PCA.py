#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[ ]:


import numpy as np
import seaborn as sns
import pandas as pd
import itertools
import copy
from scipy import stats
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import metrics
from sklearn.decomposition import PCA
from adjustText import adjust_text
import freesasa
from os import path
import os
import logomaker
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

# ## Heatmap correlation

# In[ ]:


def plot_correlation(
    self, output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Generate a correlation of each amino acid.

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

    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # calculate correlation heatmap
    dataset = _calculate_correlation(
        self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids']
    )

    # declare figure and subplots
    coeff = len(dataset.columns) / 19 * 1.05
    fig = plt.figure(figsize=(2.5 * coeff, 2.5))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax = plt.subplot(gs[0])

    # main heatmap
    heatmap = ax.pcolor(
        dataset.corr(),
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap='Greys',
        edgecolors='k',
        linewidths=0.2,
        color='darkgrey'
    )

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

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
    ax.set_xticklabels(
        list(dataset.columns),
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False
    )
    ax.set_yticklabels(
        temp_kwargs['neworder_aminoacids'],
        fontsize=6.5,
        fontname="Arial",
        color='k',
        minor=False
    )

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # _____________________________________________________________________________

    # for color bar format
    cb = plt.colorbar(
        heatmap,
        fraction=0.025,
        pad=0.05,
        aspect=5,
        ticks=[
            temp_kwargs['colorbar_scale'][0], temp_kwargs['colorbar_scale'][1]
        ],
        orientation='vertical'
    )
    cb.ax.set_yticklabels(
        cb.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k'
    )
    cb.update_ticks()
    plt.text(
        len(dataset.columns) + 1.2 * coeff,
        len(dataset.columns) / 2.5,
        'R',
        horizontalalignment='center',
        fontsize=7,
        fontname="Arial",
        color='k'
    )

    # for putting title on graph
    plt.title(
        temp_kwargs['title'],
        horizontalalignment='center',
        fontname="Arial",
        fontsize=10,
        pad=10
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()


# ## Individual correlation

# In[ ]:


def plot_individual_correlation(
    self, output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Generates a bar plot of the correlation of each amino acid mutational
    profile (row of the heatmap) with the rest of amino acids (rows)

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
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0, 1))

    # Get data
    if '*' in temp_kwargs['neworder_aminoacids']:
        temp_kwargs['neworder_aminoacids'].remove('*')
    df = _calculate_correlation(
        self.dataframe, temp_kwargs['neworder_aminoacids']
    ).mean()**2

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ticks = np.arange(0, len(df))  # label locations
    width = 0.5
    labels = temp_kwargs['neworder_aminoacids']
    # Plot figure
    ax.bar(
        ticks,
        df,
        width,
        color='blue',
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
        r'$R^2$',
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

    # show plt figure
    if temp_kwargs['show']:
        plt.show()


# ## Group correlation

# In[ ]:


def plot_group_correlation(
    self,
    r2,
    groups=['DEHKR', 'QN', 'CASTG', 'ILMV', 'WYF'],
    output=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    Determines which amino acids better represent the heatmap. Requires
    logomaker package.

    Parameters
    -----------
    self : object from class *Screen*

    r2 : float
        cutoff of the r**2 correlation value. Only values above that will be plot at the sequence logo

    groups : list, default ['DEHKR','QN','CASTG','ILMV','WYF']
        groups of aa to combine together

    output : boolean, default False

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
            
    Returns
    --------
    Use logomaker to plot the most frequent residues.
    Optional gives back the different combinations of groups and the R**2 values

    """

    # Update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # Apply parameters
    code_kwargs._parameters()

    # If there is a stop codon, delete it
    if '*' in temp_kwargs['neworder_aminoacids']:
        temp_kwargs['neworder_aminoacids'].remove('*')

    # Get R2 of each combination of amino acid substitutions
    df = _calculate_substitution_correlations(
        self, temp_kwargs['neworder_aminoacids'], groups
    )

    # Filter according the the R2 correlation value
    filtered = df.loc[df['R2'] > r2]
    logoplot = logomaker.alignment_to_matrix(list(filtered['Combinations']))

    # create Logo object
    fig = logomaker.Logo(
        logoplot,
        font_name='Arial',
        color_scheme='chemistry',
        vpad=.1,
        width=.8,
        figsize=((len(logoplot) + 1) / 2.5, 1)
    )

    # style using Logo methods
    fig.style_xticks(anchor=0, spacing=1, rotation=0)

    # No yticks and no xticks (but keep labels)
    plt.yticks([], [])
    fig.ax.tick_params(axis='both', which='both', length=0)

    # style using Axes methods
    fig.ax.set_ylabel('Bits')
    fig.ax.set_xlim([-0.5, len(logoplot) - 0.5])

    # for putting title on graph
    plt.title(
        temp_kwargs['title'],
        horizontalalignment='center',
        fontname="Arial",
        fontsize=10,
        pad=10
    )

    # save file, cannot save logo file for now

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()

    if output:
        return df


def _calculate_substitution_correlations(self, aminoacids, groups):
    '''if a set of residues was chosen, how well would they represent the entire population'''

    # Get correlation values
    corr_values = _calculate_correlation(self.dataframe, aminoacids)**2
    corr_values.reset_index(inplace=True)

    # Get combinations
    replacement_combinations = list(itertools.product(*groups))

    # Retrieve Correlation values
    df = pd.DataFrame()
    df['Aminoacids'] = list(itertools.chain.from_iterable(groups))
    for combination in replacement_combinations:  # Iterate over a combination
        temp_list = []

        # Iterate over a group of the combination
        for group, aa_selected in zip(groups, combination):
            for aa_nonselected in group:  # Find correlation values from correlation plot
                if aa_nonselected == aa_selected:
                    temp_list.append(1)
                else:
                    temp_list.append(
                        _find_correlation(
                            aa_selected, aa_nonselected, corr_values
                        )
                    )
        df[combination] = temp_list  # Store in df
    return _polishdf(df)


def _polishdf(df):
    df_mean = df.copy()
    df_mean = df.mean().to_frame()
    df_mean.reset_index(drop=False, inplace=True)
    df_mean.rename(columns={0: 'R2'}, inplace=True)
    df_mean['Combinations'] = list(df_mean['index'].apply(lambda x: ''.join(x)))
    df_mean.drop(columns=['index'], inplace=True)
    return df_mean


def _find_correlation(aa1, aa2, corr_values):
    return float(corr_values[aa1].loc[corr_values['Aminoacid'] == aa2])


# ## PCA

# In[ ]:


def plot_pca(
    self,
    mode='aminoacid',
    dimensions=[0, 1],
    adjustlabels=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    Generates a plot of two PCA dimensions.

    Parameters
    -----------
    self : object from class *Screen*

    mode : list, default 'aminoacid'
        Can also do PCA by secondary structure element if set to "secondary" or
        by individual residue if set to "individual".

    dimensions : list, default [0,1]
        Specify which two PCA dimensions to plot. By default PCA1 vs PCA2.
        Max dimension is 5.

    adjustlabels : boolean, default False
        If set to true, it will adjust the text labels so there is no overlap. It is convenient to increase
        the size of the figure, otherwise the algorithm will not find a solution. Requires to install adjustText package.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
        
        random_state : int, default 554
        
    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    """

    # load parameters
    code_kwargs._parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))

    # calculate correlation heatmap. Choose mode
    dataset = self.dataframe.copy()
    if mode == 'aminoacid':
        if '*' in temp_kwargs['neworder_aminoacids']:
            temp_kwargs['neworder_aminoacids'].remove('*')
        dataset = _calculate_correlation(
            dataset, temp_kwargs['neworder_aminoacids']
        )
        textlabels = temp_kwargs['neworder_aminoacids']
    elif mode == 'secondary':
        dataset = _calculate_correlation_bysecondary(
            dataset, self.secondary_dup
        )
        textlabels = list(dataset.columns)
    elif mode == 'individual':
        dataset = _calculate_correlation_byresidue(dataset)
        textlabels = list(dataset.columns)

    # plot using plot_clusters
    dimensionstoplot, variance = _calculate_clusters(
        dataset, dimensions, temp_kwargs['random_state']
    )

    # x and y
    x = dimensionstoplot.iloc[:, 0]
    y = dimensionstoplot.iloc[:, 1]

    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ax.scatter(x, y, s=4, c='k')

    # labels
    plt.xlabel(
        'PCA ' + str(dimensions[0] + 1) + ': ' +
        str(int(variance[dimensions[0]] * 100)) + '%',
        fontsize=10,
        labelpad=5,
        fontweight='normal'
    )
    plt.ylabel(
        'PCA ' + str(dimensions[1] + 1) + ': ' +
        str(int(variance[dimensions[1]] * 100)) + '%',
        fontsize=10,
        labelpad=-2,
        fontweight='normal'
    )

    # label of data points
    texts = _auto_text(x, y, textlabels)
    if adjustlabels is True:
        adjust_text(texts, autoalign='xy')

    # set title
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

    # show plt figure
    if temp_kwargs['show']:
        plt.show()


def _auto_text(x, y, textlabels):
    '''auto anotates text labels'''
    texts = [
        plt.annotate(
            textlabels[i],  # this is the text
            (x[i], y[i]),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(2, 2),  # distance from text to points (x,y)
            fontsize=8,
            ha='center'
        )  # horizontal alignment can be left, right or center
        for i in range(len(textlabels))
    ]
    return texts


def _calculate_clusters(dataset, dimensions, random_state):
    '''input the dataframe that needs to be correlated, the dimensions, and will calculate PCA descomposition. '''

    # call pca model
    pca = PCA(n_components=6, random_state=random_state)

    # fit model to df. use aux function correlation_aminoacids
    model = pca.fit(dataset)

    # create df with PCA data
    df_aa = pd.DataFrame((model.components_).T,
                         columns=[
                             'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6'
                         ])

    # use kmeans to cluster the two dimensions and color
    dimensionstoplot = df_aa.iloc[:, np.r_[dimensions[0], dimensions[1]]]

    return dimensionstoplot, pca.explained_variance_ratio_


def _grouby_secondary(df, secondary):
    '''
    Groups each secondary motif and makes the mean.

    Returns dataframe. Returns copy
    '''
    df = df.copy()
    df.insert(4, 'Secondary', secondary)
    df = df.groupby(['Secondary', 'Aminoacid'], as_index=False).mean()
    df = df.loc[df['Secondary'].str.startswith(('β', 'α'))]
    return df


def _calculate_correlation_bysecondary(df, secondary):
    dataset = _grouby_secondary(df, secondary)
    dataset = dataset.pivot_table(
        values='Score', index='Secondary', columns='Aminoacid'
    )
    dataset = dataset.T.corr()

    return dataset


def _calculate_correlation(df, order_aminoacids):

    dataset = df.copy()
    dataset = dataset.pivot_table(
        values='Score', index='Position', columns='Aminoacid'
    )
    dataset = dataset.corr()
    dataset = dataset.reindex(index=order_aminoacids)[order_aminoacids]

    return dataset


def _calculate_correlation_byresidue(df):

    dataset = df.copy()
    dataset = dataset.pivot_table(
        values='Score', index='Position', columns='Aminoacid'
    )
    dataset = dataset.T.corr()

    return dataset

