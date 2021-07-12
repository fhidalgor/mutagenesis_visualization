"""

"""
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
import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
import mutagenesis_visualization.main.scripts.code_utils as code_utils
from mutagenesis_visualization.main.scripts.code_heatmaps import _labels


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
        dataset = _calculate_correlation(dataset, temp_kwargs['neworder_aminoacids'])
        textlabels = temp_kwargs['neworder_aminoacids']
    elif mode == 'secondary':
        dataset = _calculate_correlation_bysecondary(dataset, self.secondary_dup)
        textlabels = list(dataset.columns)
    elif mode == 'individual':
        dataset = _calculate_correlation_byresidue(dataset)
        textlabels = list(dataset.columns)

    # plot using plot_clusters
    dimensionstoplot, variance = _calculate_clusters(dataset, dimensions, temp_kwargs['random_state'])

    # x and y
    x = dimensionstoplot.iloc[:, 0]
    y = dimensionstoplot.iloc[:, 1]

    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ax.scatter(x, y, s=4, c='k')

    # labels
    plt.xlabel(
        'PCA ' + str(dimensions[0] + 1) + ': ' + str(int(variance[dimensions[0]] * 100)) + '%',
        fontsize=10,
        labelpad=5,
        fontweight='normal'
    )
    plt.ylabel(
        'PCA ' + str(dimensions[1] + 1) + ': ' + str(int(variance[dimensions[1]] * 100)) + '%',
        fontsize=10,
        labelpad=-2,
        fontweight='normal'
    )

    # label of data points
    texts = _auto_text(x, y, textlabels)
    if adjustlabels is True:
        adjust_text(texts, autoalign='xy')

    # set title
    plt.title(temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=10, pad=5)

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
    df_aa = pd.DataFrame((model.components_).T, columns=['PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6'])

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
    dataset = dataset.pivot_table(values='Score', index='Secondary', columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset


def _calculate_correlation(df, order_aminoacids):

    dataset = df.copy()
    dataset = dataset.pivot_table(values='Score', index='Position', columns='Aminoacid')
    dataset = dataset.corr()
    dataset = dataset.reindex(index=order_aminoacids)[order_aminoacids]

    return dataset


def _calculate_correlation_byresidue(df):

    dataset = df.copy()
    dataset = dataset.pivot_table(values='Score', index='Position', columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset
