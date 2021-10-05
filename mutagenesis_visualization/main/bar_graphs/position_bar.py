"""
This module contains the class that plots the mean enrichment bar plot.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from pandas import Categorical

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import color_data


class PositionBar(Pyplot):
    """
    Class to generate a mean enrichment bar plot.
    """
    def __call__(
        self,
        position: int,
        mask_selfsubstitutions: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Choose a position and plot in a bargraph the enrichment score for each
        substitution. Red for gain of function, blue for loss of function.

        Parameters
        ----------
        position : int
            number of residue of the protein to display.

        mask_selfsubstitutions: bool, default False
            If set to true, will assing a score of 0 to each self-substitution.
            ie (A2A = 0)

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of
            the file.  Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            neworder_aminoacids: list, default list('DEKHRGNQASTPCVYMILFW*')
                Set the order (left to right) of the amino acids.

            color_gof : str, default 'red'
                Color to color mutations > 0.

            color_lof : str, default 'blue'
                Color to color mutations < 0.
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Mask self-substitutions
        self.df_output: DataFrame = self.dataframes.df_notstopcodons[replicate].copy()
        if mask_selfsubstitutions:
            self.df_output.loc[self.df_output["Sequence"] == self.df_output["Aminoacid"],"Score"] = 0

        # Select position
        self.df_output = self.df_output.loc[self.df_output['Position'] == position].copy()

        # Sort by amino acid
        self.df_output['Aminoacid'] = Categorical(self.df_output['Aminoacid'], temp_kwargs['neworder_aminoacids'])
        self.df_output = self.df_output.sort_values(by=['Aminoacid'])

        # Color
        self.df_output['Color'] = self.df_output.apply(
            color_data, axis=1, args=(
                temp_kwargs['color_gof'],
                temp_kwargs['color_lof'],
            )
        )

        # Make figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        width = 0.5

        # Color based on values
        self.ax_object.bar(
            self.df_output['Aminoacid'],
            self.df_output['Score'],
            width,
            color=self.df_output['Color'],
            ec='k',
        )

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2))
        temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')
        temp_kwargs['title_fontsize'] = kwargs.get('title_fontsize', 12)
        temp_kwargs["y_label_fontsize"] = kwargs.get('y_label_fontsize', 10)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 10)
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes parameters
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        self.ax_object.set_ylabel(
            temp_kwargs['y_label'],
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=10,
            rotation=0,
        )

        self.ax_object.set_xlabel(
            'Residue',
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            labelpad=4,
        )
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            color='k'
        )
