"""
This module contains the object to create a heatmap specifying the
selected rows.
"""
from typing import Any, Dict, Union, List, Optional
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.heatmap_utils import (
    labels,
)
from mutagenesis_visualization.main.utils.snv import add_snv_boolean
from mutagenesis_visualization.main.utils.pandas_functions import select_aa


class HeatmapRows(Pyplot):
    """
    This class plots a heatmat withe the enrichment scores.
    """
    def __call__(
        self,
        selection: Optional[List[str]] = None,
        nancolor: str = 'lime',
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a heatmap plot enrichment scores of selected aminoacids.

        Parameters
        ----------
        self : object from class *Screen*

        selection : list of aa to show (1-letter code).
            If you only want the mean displayed, type selection = 'mean'.
            By default, ["D", "K", "A", "P", "L", "W"] are displayed.

        nancolor : str, default 'lime'
            Will color np.nan values with the specified color.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        if not selection:
            selection = ["D", "K", "A", "P", "L", "W"]

        # Check if mean or not. Add group and pivot df_main.'
        if selection == 'mean':
            df_main: DataFrame = add_snv_boolean(
                self.dataframes.df_notstopcodons[replicate].copy()
            ).groupby(by='Position').mean()['Score_NaN']
            y_labels = [""]
            dataset = np.array([df_main.to_numpy()])
        else:
            df_main = select_aa(
                self.dataframes.df_stopcodons[replicate].copy(),
                selection,
                values='Score_NaN',
            )
            y_labels = list(df_main.T.columns)
            dataset = df_main.to_numpy()

        # The size can be changed. I found it empirically
        figwidth = 14 * len(dataset[0]) / 165
        figheight = 2 / 21 * len(selection)
        temp_kwargs['figsize_x'] = kwargs.get('figsize_x', figwidth)
        temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)

        self.fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
        self.gs_object: gridspec.GridSpec = gridspec.GridSpec(
            nrows=1, ncols=2, width_ratios=[len(dataset[0]), 1]
        )
        self.ax_object = plt.subplot(self.gs_object[0, 0])
        cbar1 = plt.subplot(self.gs_object[0, 1])

        # Change color of values that are NaN
        cmap = temp_kwargs['colormap']
        cmap.set_bad(color=nancolor)

        # main heatmap
        heatmap = self.ax_object.pcolormesh(
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
        self.ax_object.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
        self.ax_object.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

        # position of axis labels
        self.ax_object.tick_params('x', direction='out', pad=-2.5)
        self.ax_object.tick_params('y', direction='out', pad=0.4)

        # second axis
        ax2_object = self.ax_object.twiny()
        ax2_object.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
        ax2_object.tick_params(direction='out', pad=4)

        # Set the limits of the new axis from the original axis limits
        ax2_object.set_xlim(self.ax_object.get_xlim())

        # want a more natural, table-like display
        self.ax_object.invert_yaxis()
        self.ax_object.xaxis.tick_top()

        # so labels of x and y do not show up and my labels show up instead
        self.ax_object.set_xticklabels(
            list(self.sequence),
            fontsize=6.5,
            fontname="Arial",
            color='k',
            minor=False,
        )
        self.ax_object.set_yticklabels(
            y_labels, fontsize=6, fontname="Arial", color='k', minor=False
        )
        ax2_object.set_xticklabels(
            temp_kwargs['number_sequencelabels'][0 : len(dataset[0])],
            fontsize=10,
            fontname="Arial",
            color='k',
            minor=False
        )

        # align the labels of the y axis
        for ylabel in self.ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')

        # for coloring the residues that are 10,20...
        for xtick, color in zip(self.ax_object.get_xticklabels(),
                                temp_kwargs['color_sequencelabels']):
            xtick.set_color(color)

        # for putting title on graph
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontname="Arial",
            fontsize=temp_kwargs["title_fontsize"]
        )

        # for color bar format
        cbar1.axis('off')
        cb_object = plt.colorbar(
            heatmap,
            fraction=1,
            pad=0,
            ax=[cbar1],
            aspect=5,
            ticks=[
                temp_kwargs['colorbar_scale'][0],
                np.mean(temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]
            ],
            orientation='vertical'
        )
        cb_object.ax.set_yticklabels(
            cb_object.ax.get_yticklabels(),
            fontsize=8,
            fontname="Arial",
            color='k',
        )
        cb_object.update_ticks()
        self.gs_object.update(hspace=0.1, wspace=0.1 / len(dataset[0]) * 50)

        # remove ticks
        self.ax_object.xaxis.set_ticks_position('none')
        self.ax_object.yaxis.set_ticks_position('none')
        ax2_object.yaxis.set_ticks_position('none')
        ax2_object.xaxis.set_ticks_position('none')

        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        # load labels
        temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
        temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]
        return temp_kwargs
