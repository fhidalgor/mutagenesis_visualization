"""
This module contains the class that plots histograms.
"""
from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mutagenesis_visualization.main.classes.base_model import Pyplot


class Histogram(Pyplot):
    """
    Class to generate a kernel density plot.
    """
    def __call__(
        self,
        population: str = 'All',
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Generate a histogram plot. Can plot single nucleotide variants
        (SNVs) or non-SNVs only.

        Parameters
        ----------
        population : str, default 'All'.
            Other options are 'SNV' and 'nonSNV'.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            return_plot_object : boolean, default False
                If true, will return plotting objects (ie. fig, ax).
            bins : int, default 50.
                Number of bins for the histogram.

        Returns
        ----------
        fig : matplotlib figure and subplots
            Needs to have return_plot_object==True. By default they do
            not get returned.

        """
        temp_kwargs = self._update_kwargs(kwargs)
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])
        self.graph_parameters()

        # Select case input data
        self.dataset = self.dataframe['Score_NaN']
        if population == 'SNV':
            self.dataset = self.dataframe_snv['Score_NaN']
        elif population == 'nonSNV':
            self.dataset = self.dataframe_nonsnv['Score_NaN']

        # plot histogram
        self.ax_object = plt.hist(self.dataset, density=True, bins=temp_kwargs['bins'], color='k')

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
        temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
        temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes labels and title
        plt.xlabel(
            r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
            fontsize=temp_kwargs["x_label_fontsize"],
            fontname="Arial",
            color='k',
            labelpad=0
        )
        plt.ylabel(
            'Probability density',
            fontsize=temp_kwargs["y_label_fontsize"],
            fontname="Arial",
            color='k',
            labelpad=3
        )
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            fontname='Arial',
            color='k'
        )

        # axes limits. spacer will be 1 or the
        plt.xlim(temp_kwargs['xscale'])
        plt.xticks(
            np.arange(
                temp_kwargs['xscale'][0], temp_kwargs['xscale'][1] + temp_kwargs['tick_spacing'],
                temp_kwargs['tick_spacing']
            )
        )
        plt.ylim(temp_kwargs['yscale'])
        plt.grid()
