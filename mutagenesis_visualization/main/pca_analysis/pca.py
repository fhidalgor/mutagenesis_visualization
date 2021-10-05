"""
This module contains the pca class.
"""
from typing import Tuple, Union, Dict, Any, Literal
from pathlib import Path
import matplotlib.pyplot as plt
from adjustText import adjust_text

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pca_utils import (
    calculate_correlation,
    calculate_correlation_by_residue,
    calculate_correlation_by_secondary,
    calculate_clusters,
    auto_text,
)

MODE = Literal['aminoacid', 'secondary', 'residue']  # pylint: disable=invalid-name


class PCA(Pyplot):
    """
    This class will conduct a PCA from the enrichment scores.
    """
    def __call__(
        self,
        mode: MODE = 'aminoacid',
        dimensions: Tuple[int, int] = (0, 1),
        adjust_labels: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates a plot of two PCA dimensions.

        Parameters
        -----------
        mode : list, default 'aminoacid'
            Can also do PCA by secondary structure element if set to
            "secondary" or by individual residue if set to "residue".

        dimensions : tuple, default (0,1)
            Specify which two PCA dimensions to plot. By default PCA1 vs
            PCA2. Max dimension is 5.

        adjust_labels : boolean, default False
            If set to true, it will adjust the text labels so there is no
            overlap. It is convenient to increase the size of the figure,
            otherwise the algorithm will not find a solution. Requires to
            install adjustText package.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
            random_state : int, default 554
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # calculate correlation heatmap. Choose mode
        df_input = self.dataframes.df_notstopcodons[replicate].copy()
        if mode.lower() == 'aminoacid':
            if '*' in temp_kwargs['neworder_aminoacids']:
                temp_kwargs['neworder_aminoacids'].remove('*')
            self.df_output = calculate_correlation(df_input, temp_kwargs['neworder_aminoacids'])
            textlabels = temp_kwargs['neworder_aminoacids']
        elif mode.lower() == 'secondary':
            assert self.secondary_dup is not None, "add secondary structure information."
            self.df_output = calculate_correlation_by_secondary(df_input, self.secondary_dup)
            textlabels = list(self.df_output.columns)
        elif mode.lower() == 'individual':
            self.df_output = calculate_correlation_by_residue(df_input)
            textlabels = list(df_input.columns)

        # plot using plot_clusters
        dimensions_to_plot, variance = calculate_clusters(
            self.df_output, dimensions, temp_kwargs['random_state']
        )

        # x and y
        x = dimensions_to_plot.iloc[:, 0]
        y = dimensions_to_plot.iloc[:, 1]

        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        self.ax_object.scatter(x, y, s=4, c='k')

        # labels
        plt.xlabel(
            'PCA ' + str(dimensions[0] + 1) + ': ' + str(int(variance[dimensions[0]] * 100)) + '%',
            fontsize=temp_kwargs["x_label_fontsize"],
            labelpad=5,
            fontweight='normal'
        )
        plt.ylabel(
            'PCA ' + str(dimensions[1] + 1) + ': ' + str(int(variance[dimensions[1]] * 100)) + '%',
            fontsize=temp_kwargs["y_label_fontsize"],
            labelpad=-2,
            fontweight='normal'
        )

        # label of data points
        texts = auto_text(x, y, textlabels)
        if adjust_labels is True:
            adjust_text(texts, autoalign='xy')

        # set title
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontsize=temp_kwargs['title_fontsize'],
            pad=5
        )

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
        temp_kwargs['title_fontsize'] = kwargs.get('title_fontsize', 10)
        return temp_kwargs
