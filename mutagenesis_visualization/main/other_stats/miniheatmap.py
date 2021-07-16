"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import gridspec
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.heatmap_utils import labels
from mutagenesis_visualization.main.utils.other_stats_utils import (condense_heatmap)


class Miniheatmap(Pyplot):
    """
    Class to generate a ROC analysis.
    """
    def plot(
        self,
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Generate a miniheatmap plot enrichment scores of mutagenesis selection
        assays.

        Parameters
        ----------
        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # calculate condensed heatmap
        self.df_output = condense_heatmap(self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids'])

        # declare figure and subplots
        coeff = len(self.df_output.columns) / 19 * 1.05
        self.fig = plt.figure(figsize=(2.5 * coeff, 2.5))
        gs_object = gridspec.GridSpec(nrows=1, ncols=1)
        self.ax_object = plt.subplot(gs_object[0])

        # main heatmap
        heatmap = self.ax_object.pcolor(
            self.df_output.to_numpy(),
            vmin=temp_kwargs['colorbar_scale'][0],
            vmax=temp_kwargs['colorbar_scale'][1],
            cmap=temp_kwargs['colormap'],
            edgecolors='k',
            linewidths=0.2,
            color='darkgrey'
        )

        # for color bar format
        self.cb_object = plt.colorbar(
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
        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
       # ____________axes manipulation____________________________________________
        # put the major ticks at the middle of each cell
        self.ax_object.set_xticks(np.arange(self.df_output.shape[1]) + 0.5, minor=False)
        self.ax_object.set_yticks(np.arange(self.df_output.shape[0]) + 0.5, minor=False)

        # position of axis labels
        self.ax_object.tick_params('x', direction='out', pad=-2.5)
        self.ax_object.tick_params('y', direction='out', pad=0.4)

        # want a more natural, table-like display
        self.ax_object.invert_yaxis()
        self.ax_object.xaxis.tick_top()

        # remove ticks
        self.ax_object.xaxis.set_ticks_position('none')
        self.ax_object.yaxis.set_ticks_position('none')

        # so labels of x and y do not show up and my labels show up instead
        self.ax_object.set_xticklabels(list(self.df_output.columns), fontsize=6.5, fontname="Arial", color='k', minor=False)
        self.ax_object.set_yticklabels(
            temp_kwargs['neworder_aminoacids'], fontsize=6.5, fontname="Arial", color='k', minor=False
        )

        # align the labels of the y axis
        for ylabel in self.ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')

        # _____________________________________________________________________________

        self.cb_object.ax.set_yticklabels(self.cb_object.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k')
        self.cb_object.update_ticks()
        plt.text(
            len(self.df_output.columns) + 2,
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
