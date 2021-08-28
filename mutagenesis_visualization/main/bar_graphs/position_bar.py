"""
This module contains the class that plots the mean enrichment bar plot.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import color_data


class PositionBar(Pyplot):
    """
    Class to generate a mean enrichment bar plot.
    """
    def __call__(
        self,
        position: int,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Choose a position and plot in a bargraph the enrichment score for each
        substitution. Red for gain of function, blue for loss of function.

        Parameters
        ----------
        self : object from class *Screen*

        position : int
            number of residue of the protein to display.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of
            the file.  Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()
        # Select position
        df_output = self.dataframes.df_notstopcodons[replicate].loc[
            self.dataframes.df_notstopcodons[replicate]['Position'] == position].copy()

        # Color
        df_output['Color'] = df_output.apply(
            color_data, axis=1, args=(
                temp_kwargs['color_gof'],
                temp_kwargs['color_lof'],
            )
        )

        # make figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        width = 0.5

        # Color based on values
        self.ax_object.bar(
            df_output['Aminoacid'],
            df_output['Score'],
            width,
            color=df_output['Color'],
            ec='k',
        )

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
        temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
        temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes parameters
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        self.ax_object.set_ylabel(
            temp_kwargs['y_label'],
            fontsize=10,
            fontname="Arial",
            color='k',
            labelpad=10,
            rotation=0,
        )

        self.ax_object.set_xlabel(
            'Residue',
            fontsize=10,
            fontname="Arial",
            color='k',
            labelpad=4,
        )
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            fontname='Arial',
            color='k'
        )
