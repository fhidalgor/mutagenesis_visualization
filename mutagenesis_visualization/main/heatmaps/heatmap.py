"""
This module contains the class to plot a regular heatmap from enrichment
scores.
"""
from typing import Any, Dict, Union, Optional
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.heatmap_utils import (
    generate_cartoon,
    labels,
    hierarchical_sort,
    add_border_self_substitution,
)
from mutagenesis_visualization.main.utils.snv import add_snv_boolean
from mutagenesis_visualization.main.utils.pandas_functions import df_rearrange

EDGECOLORS: str = 'k'
LINEWIDTHS: str = 0.2
GRIDCOLOR: str = 'darkgrey'
LABELS_COLOR: str = 'k'


class Heatmap(Pyplot):
    """
    This class plots a heatmat with the enrichment scores.
    """
    def __call__(
        self,
        nancolor: str = 'lime',
        mask_selfsubstitutions: bool = False,
        color_selfsubstitutions: Optional[str] = "k",
        show_cartoon: bool = False,
        show_snv: bool = False,
        hierarchical: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a heatmap plot of the enrichment scores.

        Parameters
        ----------
        nancolor : str, default 'lime'
            Will color np.nan values with the specified color.

        mask_selfsubstitutions: bool, default False
            If set to true, will assing a score of 0 to each self-substitution.
            ie (A2A = 0)

        color_selfsubstitutions: str, default black
            If set to a color, it will color the self-substitution borders.
            Set to None to not color the self substitutions.

        show_carton : boolean, default False
            If true, the plot will display a cartoon with the secondary
            structure. The user must have added the secondary structure
            to the object.

        show_snv : boolean, default False
            If true, it will only display mutants that are a single nucleotide
            variant (SNV) of the wild-type protein sequence. The algorithm
            does not take into account the wild-type DNA allele, so it will
            include any possible mutant that is one base away.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            neworder_aminoacids: list, default list('DEKHRGNQASTPCVYMILFW*')
                Order of amino acids (y-axis) to display in the heatmap.
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # mask self-substitutions
        self.df_output: DataFrame = self.dataframes.df_stopcodons[replicate].copy()
        if mask_selfsubstitutions:
            self.df_output.loc[self.df_output["Sequence"] == self.df_output["Aminoacid"],
                               "Score_NaN"] = 0

        # sort data by rows in specified order by user
        self.df_output = df_rearrange(
            add_snv_boolean(self.df_output),
            temp_kwargs['neworder_aminoacids'],
            values='Score_NaN',
            show_snv=show_snv
        )

        # average of residues by positon
        average = [
            add_snv_boolean(self.dataframes.df_notstopcodons[replicate].copy()
                            ).groupby(by='Position').mean()['Score_NaN']
        ]
        # Create new sequence that we may to change order later
        self.sequence_updated = self.sequence

        # sort data by columns if specified by user
        if hierarchical:
            sorted_columns = hierarchical_sort(self.df_output)
            # adds the initial number so index match
            sorted_columns_corrected = sorted_columns + self.start_position
            self.df_output = self.df_output[sorted_columns_corrected]
            # had to convert back to self.df_output to be able to sort properly
            temp = DataFrame(average)
            average = temp[sorted_columns_corrected]
            # Change the sequence order
            temp_df: DataFrame = DataFrame(list(self.sequence)).T[sorted_columns]
            self.sequence_updated = list(temp_df.iloc[0])
        elif len(average) == 1:
            average = [np.array(average[0])]

        # declare figure and subplots
        figwidth = 14 * len(self.df_output.columns) / 165
        temp_kwargs['figsize_x'] = kwargs.get('figsize_x', figwidth)

        # Change parameters depending on whether cartoon is on or off
        if show_cartoon:
            figheight = 2.45
            temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)
            self.fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
            gs_object = gridspec.GridSpec(
                nrows=3,
                ncols=2,
                height_ratios=[len(self.df_output), 1, 5],
                width_ratios=[len(self.df_output.columns), 1],
            )
        else:
            figheight = 2
            temp_kwargs['figsize_y'] = kwargs.get('figsize_y', figheight)
            self.fig = plt.figure(figsize=(temp_kwargs['figsize_x'], temp_kwargs['figsize_y']))
            gs_object = gridspec.GridSpec(
                nrows=2,
                ncols=2,
                height_ratios=[len(self.df_output), 1],
                width_ratios=[len(self.df_output.columns), 1],
            )

        self.ax_object = plt.subplot(gs_object[0, 0])
        self.average_residue = plt.subplot(gs_object[1, 0])
        cbar1 = plt.subplot(gs_object[0, 1])
        cbar2 = plt.subplot(gs_object[1, 1])

        # Change color of values that are NaN
        cmap = temp_kwargs['colormap']
        cmap.set_bad(color=nancolor)

        # main heatmap
        heatmap = self.ax_object.pcolormesh(
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
                self.sequence_updated,
                temp_kwargs['neworder_aminoacids'],
                color=color_selfsubstitutions,
                lw=LINEWIDTHS * 2
            )

        # average by position
        self.average_residue.pcolormesh(
            average,
            vmin=temp_kwargs['colorbar_scale'][0],
            vmax=temp_kwargs['colorbar_scale'][1],
            cmap=cmap,
            edgecolors=EDGECOLORS,
            linewidths=LINEWIDTHS,
            antialiased=True,
            color=GRIDCOLOR
        )

        # ____________axes manipulation____________________________________________
        # put the major ticks at the middle of each cell
        self.ax_object.set_xticks(np.arange(len(self.df_output.columns)) + 0.5, minor=False)
        self.ax_object.set_yticks(np.arange(len(self.df_output)) + 0.5, minor=False)

        # position of axis labels
        self.ax_object.tick_params('x', direction='out', pad=-2.5)
        self.ax_object.tick_params('y', direction='out', pad=0.4)

        # make new axes
        self.ax_object2 = self.ax_object.twiny()
        self.ax_object3 = self.ax_object.twinx()

        # tune the axes
        self.ax_object2.set_xticks(np.arange(len(self.df_output.columns)) + 0.5, minor=False)
        self.ax_object3.set_yticks(np.arange(len(self.df_output)) + 0.5, minor=False)
        self.ax_object2.tick_params(direction='out', pad=4)
        self.ax_object3.tick_params(direction='out', pad=0.4)
        self.average_residue.tick_params(direction='out', pad=-2)
        self.average_residue.set_xticks(
            np.arange(len(self.df_output.columns)) + 0.5,
            minor=False,
        )
        self.average_residue.set_yticks(np.arange(0.5) + 0.5)

        # Set the limits of the new axis from the original axis limits
        self.ax_object2.set_xlim(self.ax_object.get_xlim())
        self.ax_object3.set_ylim(self.ax_object.get_ylim())

        # want a more natural, table-like display
        self.ax_object.invert_yaxis()
        self.ax_object.xaxis.tick_top()
        self.ax_object3.invert_yaxis()

        # remove ticks
        self.ax_object.xaxis.set_ticks_position('none')
        self.ax_object.yaxis.set_ticks_position('none')
        self.ax_object2.yaxis.set_ticks_position('none')
        self.ax_object2.xaxis.set_ticks_position('none')
        self.ax_object3.yaxis.set_ticks_position('none')
        self.average_residue.xaxis.set_ticks_position('none')
        self.average_residue.yaxis.set_ticks_position('none')

        # so labels of x and y do not show up and my labels show up instead
        self.ax_object.set_xticklabels(
            list(self.sequence_updated),
            fontsize=6.5,
            color=LABELS_COLOR,
            minor=False,
        )
        self.ax_object.set_yticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=6,
            color=LABELS_COLOR,
            minor=False,
        )
        # For numbering labels, change if hierarchical sorting is true
        if not hierarchical:
            self.ax_object2.set_xticklabels(
                temp_kwargs['number_sequencelabels'][0 : len(self.df_output.columns)],
                fontsize=10,
                color=LABELS_COLOR,
                minor=False
            )
        else:
            self.ax_object2.tick_params(direction='out', pad=7)
            self.ax_object2.set_xticklabels(
                sorted_columns_corrected,
                fontsize=5,
                color=LABELS_COLOR,
                minor=False,
                rotation=90,
            )
        self.ax_object3.set_yticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=6,
            color=LABELS_COLOR,
            minor=False,
        )
        self.average_residue.set_xticklabels(
            list(self.sequence_updated),
            fontsize=6.5,
            color=LABELS_COLOR,
            minor=False,
        )
        rowaverage = ''
        self.average_residue.set_yticklabels(
            rowaverage,
            fontsize=6,
            color=LABELS_COLOR,
            minor=False,
        )

        # align the labels of the y axis
        for ylabel in self.ax_object.get_yticklabels():
            ylabel.set_horizontalalignment('center')
        for ylabel in self.ax_object3.get_yticklabels():
            ylabel.set_horizontalalignment('center')

        # for coloring the residues that are 10,20...
        for xtick, color in zip(self.ax_object.get_xticklabels(),
                                temp_kwargs['color_sequencelabels']):
            xtick.set_color(color)
        for xtick, color in zip(self.average_residue.get_xticklabels(),
                                temp_kwargs['color_sequencelabels']):
            xtick.set_color(color)
        # _____________________________________________________________________________

        # for color bar format
        cbar1.axis('off')
        cbar2.axis('off')
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
            fontsize=6,
            color=LABELS_COLOR,
        )
        cb_object.update_ticks()
        plt.text(
            1.2 + 10 / len(self.df_output.columns),
            0.7,
            r'$\langle∆E^x_i\rangle_x$',
            transform=cbar1.transAxes,
            horizontalalignment='center',
            fontsize=7,
            color=LABELS_COLOR
        )

        gs_object.update(hspace=0.1, wspace=0.1 / len(self.df_output.columns) * 50)

        # for putting title on graph
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontsize=temp_kwargs["title_fontsize"]
        )

        # Cartoon
        if show_cartoon:
            generate_cartoon(
                self.secondary,
                self.start_position,
                gs_object,
                2,
                temp_kwargs['cartoon_colors'],
                0.025,
            )

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        # load labels
        temp_kwargs['color_sequencelabels'] = labels(self.start_position)[0]
        temp_kwargs['number_sequencelabels'] = labels(self.start_position)[1]
        return temp_kwargs
