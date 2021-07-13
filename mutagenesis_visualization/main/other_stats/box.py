"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
from numpy import array
import seaborn as sns

from mutagenesis_visualization.main.classes.base_model import Pyplot

class Box(Pyplot):
    """
    Class to generate a box plot.
    """
    def plot(self, binned_x: array, y: array, output_file: Union[None, str, Path] = None, **kwargs: Dict[str, Any],):
        """
        Generates a boxplot. Data needs to be binned prior before using
        this function.

        Parameters
        -----------
        binned_x, y : arrays
            Contain the data is going to plot.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # Make figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # Plot data
        self.ax_object = sns.boxplot(binned_x, y, color='white', fliersize=2)

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.setp(self.ax_object.artists, edgecolor='k', facecolor='w')
        plt.setp(self.ax_object.lines, color='k')

        # graph parameters
        plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k', pad=8,)
        plt.ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", color='k', labelpad=0,)
        plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k',)

        # axes limits
        plt.xlim(temp_kwargs['xscale'])
        plt.ylim(temp_kwargs['yscale'])
        plt.grid()

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
        return temp_kwargs
