"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import gridspec
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.heatmap_utils import labels
from mutagenesis_visualization.main.utils.other_stats_utils import (
    condense_heatmap, transform_dataset_offset, normalize_neighbor_effect
)


class Miniheatmap(Pyplot):
    """
    Class to generate a ROC analysis.
    """
    def __call__(
        self,
        mask_selfsubstitutions: bool = False,
        position_offset: int = 0,
        background_correction: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a miniheatmap plot enrichment scores of mutagenesis selection
        assays.

        Parameters
        ----------
        mask_selfsubstitutions: bool, default False
            If set to true, will assing a score of 0 to each self-substitution.
            i.e., (A2A = 0)

        position_offset : int, default 0
            Will group columns by residues. If the offset is not 0, it will use the values
            of the n+offset to group by. For example, you may want to see what happens when
            you have a Proline in front of the mutated residue. The algorithm can report
            the difference between the calculated value and the mean score for that particular
            substitution.
            Offset of 1 means that you evaluate the effect of following residue n+1 on n.
            Offset of -1 means that you look at the previous residue (n-1 on n).

        background_correction : boolean, default False
            If offset is nonzero, whether subtract the average effect of a substitution or not.

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

        # do offset if appropriate
        self.df_output = transform_dataset_offset(
            self.datasets[replicate],
            self.dataframes.df_notstopcodons[replicate],
            self.dataframes.df_stopcodons[replicate],
            self.sequence_raw,
            self.start_position,
            position_offset,
            stopcodons=True
        )

        # mask self-substitutions
        if mask_selfsubstitutions:
            self.df_output.loc[self.df_output["Sequence"] == self.df_output["Aminoacid"],
                               "Score"] = 0

        # calculate condensed heatmap
        self.df_output = condense_heatmap(self.df_output, temp_kwargs['neworder_aminoacids'])

        # if offset is not 0
        if background_correction and position_offset != 0:
            if '*' in temp_kwargs['neworder_aminoacids']:
                temp_kwargs['neworder_aminoacids'].remove('*')
            if '*' in temp_kwargs['aminoacids']:
                temp_kwargs['aminoacids'].remove('*')
            assert len(temp_kwargs['aminoacids']) == len(temp_kwargs['neworder_aminoacids'])

            # Add position_offset sequence to df
            self.df_output = transform_dataset_offset(
                self.datasets[replicate],
                self.dataframes.df_notstopcodons[replicate],
                self.dataframes.df_stopcodons[replicate],
                self.sequence_raw,
                self.start_position,
                position_offset,
                stopcodons=False
            )
            # mask self-substitutions
            if mask_selfsubstitutions:
                self.df_output.loc[self.df_output["Sequence"] == self.df_output["Aminoacid"],
                                   "Score"] = 0

            # calculate mean effect using condensed heatmap
            df_condensed_heatmap = condense_heatmap(
                self.df_output, temp_kwargs['neworder_aminoacids']
            )
            self.df_output = normalize_neighbor_effect(
                self.df_output, df_condensed_heatmap, temp_kwargs['neworder_aminoacids'],
                temp_kwargs['aminoacids']
            )

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
        self.ax_object.set_xticklabels(
            list(self.df_output.columns), fontsize=temp_kwargs["x_label_fontsize"], color='k', minor=False
        )
        self.ax_object.set_yticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            minor=False
        )

        # align the labels of the y axis
        for ylabel in self.ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')

        # _____________________________________________________________________________

        self.cb_object.ax.set_yticklabels(
            self.cb_object.ax.get_yticklabels(), fontsize=7, color='k'
        )
        self.cb_object.update_ticks()
        plt.text(
            len(self.df_output.columns) + 2,
            7.8,
            r'$\langle∆E^x_i\rangle_x$',
            horizontalalignment='center',
            fontsize=7,
            color='k'
        )

        # for putting title on graph
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontsize=temp_kwargs["title_fontsize"],
            pad=10
        )
        plt.ylabel('Amino Acid Substitution', fontsize=temp_kwargs["y_label_fontsize"], labelpad=-1)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        # load labels
        temp_kwargs['aminoacids'] = self.aminoacids
        temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
        temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]
        temp_kwargs["y_label_fontsize"] = kwargs.get('y_label_fontsize', 6.5)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 6.5)
        return temp_kwargs
