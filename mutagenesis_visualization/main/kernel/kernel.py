"""
This module contains the kernel class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

from mutagenesis_visualization.main.classes.base_model import Pyplot


class Kernel(Pyplot):
    """
    Class to generate a kernel density plot.
    """
    def __call__(
        self,
        cumulative: bool = False,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
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
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        self.fig = plt.figure(figsize=temp_kwargs['figsize'])

        # plot kernel
        self.ax_object = sns.kdeplot(
            self.dataframe['Score_NaN'],
            cumulative=cumulative,
            color='red',
            lw=2,
        )
        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.xlabel(
            temp_kwargs['x_label'],
            fontsize=temp_kwargs["x_label_fontsize"],
            fontname='Arial',
            color='k',
            labelpad=0
        )
        plt.ylabel(
            temp_kwargs['y_label'],
            fontsize=temp_kwargs["y_label_fontsize"],
            fontname='Arial',
            color='k',
            labelpad=3
        )
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            fontname='Arial',
            color='k'
        )
        plt.xlim(temp_kwargs['xscale'])
        plt.grid()
        self.ax_object.get_legend().remove()

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
        temp_kwargs['x_label'] = kwargs.get('x_label', r'$∆E^i_x$')
        temp_kwargs['y_label'] = kwargs.get('y_label', 'Probability density')
        return temp_kwargs