"""
This module contains the correlation class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.heatmap_utils import labels
from mutagenesis_visualization.main.utils.pca_utils import calculate_correlation


class Correlation(Pyplot):
    """
    This class will conduct a correlation from the enrichment scores.
    """
    def __call__(
        self,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a correlation of each amino acid.

        Parameters
        ----------
        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            colorbar_scale: tuple, default (-1, 1)
                Scale min and max used in heatmaps and correlation heatmaps.
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # calculate correlation heatmap
        self.df_output = calculate_correlation(
            self.dataframes.df_stopcodons[replicate],
            temp_kwargs['neworder_aminoacids'],
        )

        # declare figure and subplots
        self.fig = plt.figure(figsize=(2.5 * len(self.df_output.columns) / 19 * 1.05, 2.5))
        self.gs_object = gridspec.GridSpec(nrows=1, ncols=1)
        self.ax_object = plt.subplot(self.gs_object[0])

        # main heatmap
        heatmap = self.ax_object.pcolor(
            self.df_output.corr(),
            vmin=temp_kwargs['colorbar_scale'][0],
            vmax=temp_kwargs['colorbar_scale'][1],
            cmap='Greys',
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
            ticks=[temp_kwargs['colorbar_scale'][0], temp_kwargs['colorbar_scale'][1]],
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
        self.ax_object.set_xticklabels(
            list(self.df_output.columns),
            fontsize=6.5,
            color='k',
            minor=False,
        )
        self.ax_object.set_yticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=6.5,
            color='k',
            minor=False,
        )

        # align the labels of the y axis
        for ylabel in self.ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')

        # _____________________________________________________________________________

        self.cb_object.ax.set_yticklabels(
            self.cb_object.ax.get_yticklabels(),
            fontsize=7,
            color='k',
        )
        self.cb_object.update_ticks()
        plt.text(
            len(self.df_output.columns) + 1.2 * len(self.df_output.columns) / 19 * 1.05,
            len(self.df_output.columns) / 2.5,
            'R',
            horizontalalignment='center',
            fontsize=7,
            color='k'
        )

        # for putting title on graph
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontsize=10,
            pad=10,
        )

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
        temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]
        return temp_kwargs
