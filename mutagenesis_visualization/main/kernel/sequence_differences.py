"""
This module will host the sequence_differences class
"""

from typing import Tuple, Union, Dict, Any, List, TYPE_CHECKING, Optional
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import get_fitness_changes
if TYPE_CHECKING:
    from mutagenesis_visualization.main.classes.screen import Screen


class SequenceDifferences(Pyplot):
    """
    Class to generate the sequence differences plot.
    """
    def __call__(
        self,
        screen_object: 'Screen',
        map_sequence_changes: List[Tuple[int, int]],
        legend_labels: Optional[Tuple[str, str]] = None,
        replicate: int = -1,
        replicate_second_object: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate two histogram plots. The first plot will have the impact
        on fitness to go from protein A -> B, and the second plot will
        contain the B -> A effect.

        Parameters
        ----------
        screen_object : *Screen* object or list containing *Screen*
            objects.

        map_sequence_changes: list of tuples
            Set the residues that differ between protein A and protein B.
            Example: [(1, 1), (12, 12), (15, 16)]. In the example, the
            algorithm will compare the residue 1 and 12 of each protein,
            and the residue 15 of protein A vs the residue 16 of protein B.

        legend_labels: tuple of str
            Set the labels of the legend.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        replicate_second_object : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            bins : int or str, default 'auto'.
                Number of bins for the histogram. By default it will
                automatically decide the number of bins.

        """
        temp_kwargs = self._update_kwargs(kwargs)
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])
        self.graph_parameters()
        if not legend_labels:
            legend_labels = ("A -> B", "B -> A")

        # Select case input data
        self.df_output = get_fitness_changes(
            map_sequence_changes, self.dataframes.df_notstopcodons[replicate],
            screen_object.dataframes.df_notstopcodons[replicate_second_object]
        )

        # plot histogram
        self.ax_object = plt.hist(
            self.df_output["A_to_B"],
            density=True,
            bins=temp_kwargs['bins'],
            color='blue',
            label=legend_labels[0],
            alpha=0.5
        )
        self.ax_object = plt.hist(
            self.df_output["B_to_A"],
            density=True,
            bins=temp_kwargs['bins'],
            color='red',
            label=legend_labels[1],
            alpha=0.5
        )
        plt.legend(loc='best', frameon=False, handlelength=1, handletextpad=0.5)

        # legend????
        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4, 2))
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes labels and title
        plt.xlabel(
            r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            labelpad=0
        )
        plt.ylabel(
            'Probability density',
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=3
        )
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            color='k'
        )

        # axes limits. spacer will be 1 or the
        if temp_kwargs['xscale'] != (None, None):
            plt.xlim(temp_kwargs['xscale'])
            plt.xticks(
                np.arange(
                    temp_kwargs['xscale'][0],
                    temp_kwargs['xscale'][1] + temp_kwargs['tick_spacing'],
                    temp_kwargs['tick_spacing']
                )
            )
        if temp_kwargs['yscale'] != (None, None):
            plt.ylim(temp_kwargs['yscale'])
        plt.grid()
