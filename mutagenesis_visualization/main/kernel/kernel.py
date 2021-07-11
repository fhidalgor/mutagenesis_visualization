"""
This module contains the kernel class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
import seaborn as sns

from mutagenesis_visualization.main.classes.base_model import Pyplot


class Kernel(Pyplot):
    """
    Class to generate
    """
    def plot(
        self,
        cumulative=False,
        output_file: Union[None, str, Path] = None,
        **kwargs,
    ):
        """
        Plot univariate or bivariate distributions using kernel density estimation.
        Parameters
        ----------
        cumulative : bool, optional, default False
            If True, estimate a cumulative distribution function.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of the file.
            Example: 'path/filename.png' or 'path/filename.svg'.
        **kwargs : other keyword arguments
            return_plot_object : boolean, default False
                If true, will return plotting objects (ie. fig, ax_object).
        Returns
        ----------
        fig, ax_object : matplotlib figure and subplots
            Needs to have return_plot_object=True. By default they do
            not get returned.
        """
        temp_kwargs = self._update_kwargs(kwargs)
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])
        self._load_parameters()

        # plot
        self.ax_object = sns.kdeplot(
            self.dataset,
            cumulative=cumulative,
            color='red',
            lw=2,
        )
        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _tune_plot(self, temp_kwargs) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname='Arial', color='k', labelpad=0)
        plt.ylabel(['y_label'], fontsize=10, fontname='Arial', color='k', labelpad=3)
        plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
        plt.xlim(temp_kwargs['xscale'])
        plt.grid()

    def _update_kwargs(self, kwargs) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
        temp_kwargs['x_label'] = kwargs.get('x_label', r'$∆E^i_x$')
        temp_kwargs['y_label'] = kwargs.get('y_label', 'Probability density')
        return temp_kwargs

    def return_plot_object(self,):
        """
        Return matplotlib object.
        """
        return self.fig, self.ax_object
