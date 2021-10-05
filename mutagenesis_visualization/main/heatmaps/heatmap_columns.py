"""
This module contains the object to create a heatmap specifying the
selected columns.
"""
from typing import Any, Dict, Tuple, Union, Optional
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.heatmaps.heatmap import (
    LABELS_COLOR, LINEWIDTHS, GRIDCOLOR, EDGECOLORS
)
from mutagenesis_visualization.main.utils.heatmap_utils import labels, add_border_self_substitution
from mutagenesis_visualization.main.utils.pandas_functions import df_rearrange


class HeatmapColumns(Pyplot):
    """
    This class plots a heatmap with the enrichment scores where you
    can show selected columns.
    """
    def __call__(
        self,
        segment: Tuple[int, int],
        ylabel: bool = True,
        nancolor: str = 'lime',
        mask_selfsubstitutions: bool = False,
        color_selfsubstitutions: Optional[str] = "k",
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a heatmap plot enrichment scores but only plots a selected segment.

        Parameters
        ----------
        segment : Tuple[int]
            Segment is typed as [20,40] and includes both residues 20 and 40.

        ylabel : str, default True
            Choose False to hide.

        nancolor : str, default 'lime'
            Will color np.nan values with the specified color.

        mask_selfsubstitutions: bool, default False
            If set to true, will assing a score of 0 to each self-substitution.
            ie (A2A = 0)

        color_selfsubstitutions: str, default black
            If set to a color, it will color the self-substitution borders.
            Set to None to not color the self substitutions.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of the file.
            Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # mask self-substitutions
        df_main: DataFrame = self.dataframes.df_stopcodons[replicate].copy()
        if mask_selfsubstitutions:
            df_main.loc[df_main["Sequence"] == df_main["Aminoacid"], "Score_NaN"] = 0

        # sort data in specified order by user
        df_main = df_rearrange(df_main, temp_kwargs['neworder_aminoacids'], values='Score_NaN')

        # select subset
        self.df_output = df_main.iloc[:, segment[0] - self.start_position : segment[1] -
                                      self.start_position + 1]

        # the size can be changed
        figwidth = 14 * len(self.df_output.columns) / 165
        figheight = 2
        temp_kwargs['figsize_x'] = kwargs.get('figsize_x', figwidth)
        temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)

        self.fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
        gs_object: gridspec.GridSpec = gridspec.GridSpec(nrows=1, ncols=1)
        # needed to set autoscale off to avoid missalignment
        self.ax_object = plt.subplot(gs_object[0])

        # Change color of values that are NaN
        cmap = temp_kwargs['colormap']
        cmap.set_bad(color=nancolor)

        # main heatmap
        self.ax_object.pcolormesh(
            self.df_output,
            vmin=temp_kwargs['colorbar_scale'][0],
            vmax=temp_kwargs['colorbar_scale'][1],
            cmap=cmap,
            edgecolors=EDGECOLORS,
            linewidths=LINEWIDTHS,
            antialiased=True,
            color=GRIDCOLOR
        )

        # add border to self-substitution square
        if color_selfsubstitutions:
            add_border_self_substitution(
                self.ax_object,
                self.sequence[segment[0] - self.start_position : segment[1] - self.start_position +
                              1],
                temp_kwargs['neworder_aminoacids'],
                color=color_selfsubstitutions,
                lw=LINEWIDTHS * 2
            )

        # put the major ticks at the middle of each cell
        self.ax_object.set_xticks(
            np.arange(len(self.df_output.columns)) + 0.5,
            minor=False,
        )
        self.ax_object.set_yticks(np.arange(len(self.df_output)) + 0.5, minor=False)

        # position of axis labels
        self.ax_object.tick_params('x', direction='out', pad=-2.5)
        self.ax_object.tick_params('y', direction='out', pad=0.4)

        # second axis
        ax2_object = self.ax_object.twiny()
        ax2_object.set_xticks(np.arange(len(self.df_output.columns)) + 0.5, minor=False)
        ax2_object.tick_params(direction='out', pad=4)

        # Set the limits of the new axis from the original axis limits
        ax2_object.set_xlim(self.ax_object.get_xlim())

        # want a more natural, table-like display
        self.ax_object.invert_yaxis()
        self.ax_object.xaxis.tick_top()

        # so labels of x and y do not show up and my labels show up instead
        self.ax_object.set_xticklabels(
            list(self.sequence
                 )[segment[0] - self.start_position : segment[1] - self.start_position + 1],
            fontsize=temp_kwargs["x_label_fontsize"],
            color=LABELS_COLOR,
            minor=False,
        )
        self.ax_object.set_yticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=temp_kwargs["y_label_fontsize"],
            color=LABELS_COLOR,
            minor=False,
        )
        if not ylabel:
            self.ax_object.set_yticklabels([])


        ax2_label = (segment[1] - segment[0] + 1) * ['']
        ax2_label[0] = str(segment[0])
        ax2_label[-1] = str(segment[1])
        ax2_object.set_xticklabels(ax2_label, fontsize=7, color=LABELS_COLOR, minor=False)

        # align the labels of the y axis
        for ylabel in self.ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')

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
        temp_kwargs["y_label_fontsize"] = kwargs.get('y_label_fontsize', 6)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 6.5)
        return temp_kwargs
