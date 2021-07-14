"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from matplotlib import ticker
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.other_stats_utils import (
    select_grouping,
    merge_class_variants,
    roc_auc,
)
from mutagenesis_visualization.main.utils.heatmap_utils import labels


class Miniheatmap(Pyplot):
    """
    Class to generate a ROC analysis.
    """
    def plot(
        self,
        offset: float = 0,
        background_correction: bool = True,
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ):
        """
        Generate a miniheatmap plot enrichment scores of mutagenesis selection
        assays.

        Parameters
        ----------
        self : object from class *Screen*

        offset : int, default 0
            Will group columns by residues. If the offset is not 0, it will
            use the values of the n+offset to group by. For example, you
            may want to see what happens when you have a Proline in front
            of the mutated residue. The algorithm can report the
            difference between the calculated value and the mean score
            for that particular substitution. Offset of 1 means that you
            evaluate the effect of following residue n+1 on n. Offset of
            -1 means that you look at the previous residue (n-1 on n).

        background_correction : boolean, default True
            If offset is nonzero, whether subtract the average effect of a
            substitution or not.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # do offset if appropriate
        dataframe_stopcodons = transform_dataset_offset(self, offset)

        # calculate condensed heatmap
        dataset = condense_heatmap(dataframe_stopcodons, temp_kwargs['neworder_aminoacids'])

        # if offset is not 0
        if background_correction and offset != 0:
            if '*' in temp_kwargs['neworder_aminoacids']:
                temp_kwargs['neworder_aminoacids'].remove('*')
            # do offset, no stop codons
            dataset = _normalize_neighboreffect(self, offset, temp_kwargs['neworder_aminoacids'])

        fig, ax, cb = _plot_miniheatmap(dataset, output_file, temp_kwargs)

        # return matplotlib object
        if temp_kwargs['return_plot_object']:
            return fig, ax, cb

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        # load labels
        temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
        temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]
        return temp_kwargs


def _plot_miniheatmap(df, output_file, temp_kwargs):
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
    coeff = len(df.columns) / 19 * 1.05
    fig = plt.figure(figsize=(2.5 * coeff, 2.5))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax = plt.subplot(gs[0])

    # main heatmap
    heatmap = ax.pcolor(
        df.to_numpy(),
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=temp_kwargs['colormap'],
        edgecolors='k',
        linewidths=0.2,
        color='darkgrey'
    )

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
    ax.set_xticklabels(list(df.columns), fontsize=6.5, fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(
        temp_kwargs['neworder_aminoacids'], fontsize=6.5, fontname="Arial", color='k', minor=False
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
            temp_kwargs['colorbar_scale'][0],
            np.mean(temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]
        ],
        orientation='vertical'
    )
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k')
    cb.update_ticks()
    plt.text(
        len(df.columns) + 2,
        7.8,
        r'$\langle∆E^x_i\rangle_x$',
        horizontalalignment='center',
        fontsize=7,
        fontname="Arial",
        color='k'
    )

    # for putting title on graph
    plt.title(
        temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=10, pad=10
    )
    plt.ylabel('Amino Acid Substitution', fontsize=10, labelpad=-1)

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    return fig, ax, cb
