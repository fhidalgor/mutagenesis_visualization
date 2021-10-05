"""
This module contains the class that plots the mean position bar plot.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

from mutagenesis_visualization.main.bar_graphs.mean_counts import input_text
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import color_data
from mutagenesis_visualization.main.utils.heatmap_utils import generate_cartoon


class EnrichmentBar(Pyplot):
    """
    Class to generate a enrichment bar plot per position.
    """
    def __call__(
        self,
        mode: str = "mean",
        show_cartoon: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot in a bargraph the enrichment for each residue of the
        protein. Red for gain of function, blue for loss of function.

        Parameters
        ----------
        mode : str, default 'mean'
            Specify what enrichment scores to show. If mode = 'mean', it
            will show the mean of each position. If mode = 'A', it will
            show the alanine substitution profile. Can be used for each
            amino acid. Use the one-letter code and upper case.

        show_cartoon : boolean, default False
            If true, the plot will display a cartoon with the secondary
            structure. The user must have added the secondary structure
            to the object.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            color_gof : str, default 'red'
                Choose color to color positions with an enrichment score > 0.

            color_lof : str, default 'blue'
                Choose color to color positions with an enrichment score < 0.
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Select grouping
        if mode.lower() == 'mean':
            self.df_output = self.dataframes.df_notstopcodons[replicate].groupby(
                'Position', as_index=False
            ).mean()
        else:
            self.df_output = self.dataframes.df_notstopcodons[replicate].loc[
                self.dataframes.df_notstopcodons[replicate]['Aminoacid'] == mode].copy()

        self.df_output['Color'] = self.df_output.apply(
            color_data, axis=1, args=(
                temp_kwargs['color_gof'],
                temp_kwargs['color_lof'],
            )
        )

        # make figure
        if show_cartoon:
            self.fig = plt.figure(figsize=temp_kwargs['figsize'])
            gs_object: gridspec.GridSpec = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
            self.ax_object = plt.subplot(gs_object[0])
        else:
            self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # Color based on values
        self.ax_object.bar(
            self.df_output['Position'],
            self.df_output['Score'],
            width=1.2,
            color=self.df_output['Color'],
            snap=False,
        )

        # cartoon
        title_pad: float = 0
        if show_cartoon:
            generate_cartoon(
                self.secondary,
                self.start_position,
                gs_object,
                1,
                temp_kwargs['cartoon_colors'],
                bottom_space=-0.78,
                show_labels=False,
            )
            title_pad = 2.5

        # Plot title
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            color='k',
            pad=title_pad
        )

        # Put text labels
        input_text(temp_kwargs['text_labels'])

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (15, 2.5))
        temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')
        temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 10)
        temp_kwargs["y_label_fontsize"] = kwargs.get('y_label_fontsize', 10)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 10)
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # add grid
        if temp_kwargs['grid']:
            self.ax_object.grid(which='major', color='silver', linewidth=1)
            self.ax_object.grid(which='minor', color='silver', linewidth=0.25)
            # Show the minor ticks and grid.
            self.ax_object.minorticks_on()

        # axes parameters
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        self.ax_object.set_ylabel(
            temp_kwargs['y_label'],
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=10,
            rotation=0
        )
        self.ax_object.set_xticks(
            np.arange(
                self.start_position,
                len(self.df_output) + self.start_position,
                temp_kwargs['tick_spacing'],
            )
        )
        self.ax_object.set_xlabel('Residue', fontsize=temp_kwargs["x_label_fontsize"], color='k', labelpad=4)
        self.ax_object.set_xlim(
            self.start_position - 0.1,
            len(self.df_output) + self.start_position - 1 + 0.1,
        )
