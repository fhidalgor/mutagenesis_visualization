"""
This module contains the class rank.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from mutagenesis_visualization.main.classes.base_model import Pyplot


class Rank(Pyplot):
    """
    Class to generate a mean enrichment bar plot.
    """
    def __call__(
        self,
        mode: str = 'pointmutant',
        output_file: Union[None, str, Path] = None,
        replicate: int = -1,
        **kwargs: Any
    ) -> None:
        """
        Generate a rank plot so every mutation/residue is sorted based
        on enrichment score.

        Parameters
        ----------
        mode : str, default 'pointmutant'.
            Alternative set to "mean" for the mean of each position

        outdf : boolean, default False
            If set to true, will return the df with the rank of mutations

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of the file.
            Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Sort by enrichment scores
        self.df_output = self.dataframes.df_notstopcodons[replicate].sort_values(by=['Score']
                                                                                 ).copy()

        # Chose mode:
        if mode == 'mean':
            self.df_output = self.df_output.groupby(by=['Position'], as_index=False).mean()
            self.df_output.sort_values(by=['Score'], inplace=True)

        # create figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4, 2))
        temp_kwargs['x_label'] = kwargs.get('x_label', 'Rank')
        temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')
        temp_kwargs['title_fontsize'] = kwargs.get('title_fontsize', 12)
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Scatter data points
        plt.scatter(np.arange(len(self.df_output), 0, -1), self.df_output['Score'], c='k', s=1)

        # Titles
        plt.title(temp_kwargs['title'], fontsize=temp_kwargs['title_fontsize'], color='k', pad=8)
        # Labels
        plt.ylabel(
            temp_kwargs['y_label'],
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=0
        )
        plt.xlabel(
            temp_kwargs['x_label'],
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k'
        )

        # other graph parameters
        plt.xlim(temp_kwargs['xscale'])
        plt.ylim(temp_kwargs['yscale'])
