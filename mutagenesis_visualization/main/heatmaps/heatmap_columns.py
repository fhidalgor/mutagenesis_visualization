"""
This module contains the object to create a heatmap specifying the
selected columns.
"""
from typing import Any, Dict, Tuple, Union, Optional
from pathlib import Path
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.heatmap_utils import labels
from mutagenesis_visualization.main.utils.pandas_functions import df_rearrange


class HeatmapColumns(Pyplot):
    """
    This class plots a heatmap with the enrichment scores where you
    can show selected columns.
    """

    def __call__(
        self,
        segment: Tuple[int],
        ylabel_color: str = 'k',
        nancolor: str = 'lime',
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ):
        """
        Generate a heatmap plot enrichment scores but only plots a selected segment.

        Parameters
        ----------
        self : object from class *Screen*

        segment : Tuple[int]
            Segment is typed as [20,40] and includes both residues 20 and 40.

        ylabel_color : str, default 'k'
            Choose white if you don't want amino acid y axis label.

        nancolor : str, default 'lime'
            Will color np.nan values with the specified color.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of the file.
            Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # sort data in specified order by user
        df_whole: pd.DataFrame = df_rearrange(
            self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids'], values='Score_NaN'
        )

        # select subset
        self.df_output = df_whole.iloc[:, segment[0] - self.start_position : segment[1] -
                                       self.start_position + 1]

        # the size can be changed
        figwidth = 2 * len(self.df_output.columns) / 22
        figheight = 2
        temp_kwargs['figsize_x'] = kwargs.get('figsize_x', figwidth)
        temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)

        self.fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
        gs_object = gridspec.GridSpec(nrows=1, ncols=1)
        # needed to set autoscale off to avoid missalignment
        ax_object = plt.subplot(gs_object[0])

        # Change color of values that are NaN
        cmap = temp_kwargs['colormap']
        cmap.set_bad(color=nancolor)

        # main heatmap
        ax_object.pcolormesh(
            self.df_output,
            vmin=temp_kwargs['colorbar_scale'][0],
            vmax=temp_kwargs['colorbar_scale'][1],
            cmap=cmap,
            edgecolors='k',
            linewidths=0.2,
            antialiased=True,
            color='darkgrey'
        )

        # put the major ticks at the middle of each cell
        ax_object.set_xticks(
            np.arange(len(self.df_output.columns)) + 0.5,
            minor=False,
        )
        ax_object.set_yticks(np.arange(len(self.df_output)) + 0.5, minor=False)

        # position of axis labels
        ax_object.tick_params('x', direction='out', pad=-2.5)
        ax_object.tick_params('y', direction='out', pad=0.4)

        # second axis
        ax2_object = ax_object.twiny()
        ax2_object.set_xticks(np.arange(len(self.df_output.columns)) + 0.5, minor=False)
        ax2_object.tick_params(direction='out', pad=4)

        # Set the limits of the new axis from the original axis limits
        ax2_object.set_xlim(ax_object.get_xlim())

        # want a more natural, table-like display
        ax_object.invert_yaxis()
        ax_object.xaxis.tick_top()

        # so labels of x and y do not show up and my labels show up instead
        ax_object.set_xticklabels(
            list(self.sequence
                 )[segment[0] - self.start_position : segment[1] - self.start_position + 1],
            fontsize=6.5,
            fontname="Arial",
            color='k',
            minor=False,
        )
        ax_object.set_yticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=6,
            fontname="Arial",
            color=ylabel_color,
            minor=False,
        )

        ax2_label = (segment[1] - segment[0] + 1) * ['']
        ax2_label[0] = segment[0]
        ax2_label[-1] = segment[1]
        ax2_object.set_xticklabels(ax2_label, fontsize=7, fontname="Arial", color='k', minor=False)

        # align the labels of the y axis
        for ylabel in ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')

        # remove ticks
        ax_object.xaxis.set_ticks_position('none')
        ax_object.yaxis.set_ticks_position('none')
        ax2_object.yaxis.set_ticks_position('none')
        ax2_object.xaxis.set_ticks_position('none')

        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _update_kwargs(self, kwargs) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] =  super()._update_kwargs(kwargs)
        # load labels
        temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
        temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]
        return temp_kwargs
