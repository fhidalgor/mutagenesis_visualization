"""
This module contains the correlation class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pca_utils import calculate_correlation


class IndividualCorrelation(Pyplot):
    """
    This class will conduct an individual correlation from the enrichment
    scores.
    """
    def __call__(
        self,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates a bar plot of the correlation of each amino acid mutational
        profile (row of the heatmap) with the rest of amino acids (rows)

        Parameters
        -----------
        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Get data
        if '*' in temp_kwargs['neworder_aminoacids']:
            temp_kwargs['neworder_aminoacids'].remove('*')
        self.df_output = calculate_correlation(
            self.dataframes.df_notstopcodons[replicate], temp_kwargs['neworder_aminoacids']
        ).mean()**2

        # Make figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        ticks = np.arange(0, len(self.df_output))  # label locations
        width = 0.5
        # Plot figure
        self.ax_object.bar(
            ticks,
            self.df_output,
            width,
            color='blue',
            ec='k',
        )

        # graph parameters
        self.ax_object.set_xticks(ticks)
        self.ax_object.set_xticklabels(
            temp_kwargs['neworder_aminoacids'],
            fontsize=9,
            color='k',
            minor=False,
            rotation=0,
        )
        self.ax_object.set_ylabel(
            r'$R^2$',
            fontsize=10,
            color='k',
            labelpad=12,
            rotation=0,
        )
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontsize=10,
            pad=5,
        )

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
        temp_kwargs['yscale'] = kwargs.get('yscale', (0, 1))
        return temp_kwargs
