"""
This module contains the class that plots scatters.
"""
from typing import Union, Dict, Any
from pathlib import Path
import numpy as np
from pandas.core.frame import DataFrame
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import ticker
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import (
    process_mean_residue,
    process_by_pointmutant,
)


class Scatter(Pyplot):
    """
    Class to generate a kernel density plot.
    """
    def __call__(
        self,
        screen_object: Any,
        mode: str = 'pointmutant',
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Generate a scatter plot between object and a second object of the
        same class.

        Parameters
        ----------
        screen_object : object from class *Screen* to do the scatter with

        mode : str, default 'pointmutant'.
            Alternative set to "mean" for the mean of each position.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Chose mode:
        if mode.lower() == 'pointmutant':
            df_output: DataFrame = process_by_pointmutant(self.dataframe, screen_object.dataframe)
        else:
            df_output = process_mean_residue(self.dataframe, screen_object.dataframe)

        # create figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # Scatter data points
        plt.scatter(
            df_output['dataset_1'],
            df_output['dataset_2'],
            c='k',
            s=8,
            alpha=0.5,
            rasterized=True,
            label='_nolegend_'
        )

        # correlation
        _, _, r_value, _, _ = linregress(df_output['dataset_1'], df_output['dataset_2'])

        # fit and graph line
        fit = np.polyfit(df_output['dataset_1'], df_output['dataset_2'], 1)
        plt.plot(
            np.unique(df_output['dataset_1']),
            np.poly1d(fit)(np.unique(df_output['dataset_1'])),
            color='r',
            linewidth=1,
            label="$R^2$ = {}".format(str(round(r_value**2, 2)))
        )

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Titles
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            fontname='Arial',
            color='k',
            pad=8
        )
        plt.ylabel(
            temp_kwargs['y_label'],
            fontsize=temp_kwargs["y_label_fontsize"],
            fontname="Arial",
            color='k',
            labelpad=0
        )
        plt.xlabel(
            temp_kwargs['x_label'],
            fontsize=temp_kwargs["x_label_fontsize"],
            fontname="Arial",
            color='k'
        )

        plt.grid()

        # other graph parameters
        plt.xlim(temp_kwargs['xscale'])
        plt.ylim(temp_kwargs['yscale'])
        self.ax_object.xaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))
        self.ax_object.yaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))
        plt.gca().set_aspect('equal', adjustable='box')
        plt.draw()

        # Legend
        plt.legend(loc='upper left', handlelength=0, handletextpad=0, frameon=False, fontsize=10)
