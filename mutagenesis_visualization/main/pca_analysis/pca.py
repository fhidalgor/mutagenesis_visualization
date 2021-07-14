"""
This module contains the pca class.
"""
from typing import List, Union, Dict, Any
import copy
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


class PCA(Pyplot):
    """
    This class will conduct a PCA from the enrichment scores.
    """
    def plot(
        self,
        mode: str = 'aminoacid',
        dimensions: List[int] = [0, 1],
        adjust_labels: bool = False,
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ):
        """
        Generates a plot of two PCA dimensions.

        Parameters
        -----------
        self : object from class *Screen*

        mode : list, default 'aminoacid'
            Can also do PCA by secondary structure element if set to
            "secondary" or by individual residue if set to "individual".

        dimensions : list, default [0,1]
            Specify which two PCA dimensions to plot. By default PCA1 vs
            PCA2. Max dimension is 5.

        adjust_labels : boolean, default False
            If set to true, it will adjust the text labels so there is no
            overlap. It is convenient to increase the size of the figure,
            otherwise the algorithm will not find a solution. Requires to
            install adjustText package.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
            random_state : int, default 554
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # calculate correlation heatmap. Choose mode
        dataset = self.dataframe.copy()
        if mode == 'aminoacid':
            if '*' in temp_kwargs['neworder_aminoacids']:
                temp_kwargs['neworder_aminoacids'].remove('*')
            dataset = calculate_correlation(dataset, temp_kwargs['neworder_aminoacids'])
            textlabels = temp_kwargs['neworder_aminoacids']
        elif mode == 'secondary':
            dataset = calculate_correlation_by_secondary(dataset, self.secondary_dup)
            textlabels = list(dataset.columns)
        elif mode == 'individual':
            dataset = calculate_correlation_by_residue(dataset)
            textlabels = list(dataset.columns)

        # plot using plot_clusters
        dimensions_to_plot, variance = calculate_clusters(
            dataset, dimensions, temp_kwargs['random_state']
        )

        # x and y
        x = dimensions_to_plot.iloc[:, 0]
        y = dimensions_to_plot.iloc[:, 1]

        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        self.ax_object.scatter(x, y, s=4, c='k')

        # labels
        plt.xlabel(
            'PCA ' + str(dimensions[0] + 1) + ': ' + str(int(variance[dimensions[0]] * 100)) + '%',
            fontsize=10,
            labelpad=5,
            fontweight='normal'
        )
        plt.ylabel(
            'PCA ' + str(dimensions[1] + 1) + ': ' + str(int(variance[dimensions[1]] * 100)) + '%',
            fontsize=10,
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
            fontname="Arial",
            fontsize=10,
            pad=5
        )

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
        return temp_kwargs
