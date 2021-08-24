"""
This module contains the multiple kernel class.
"""

from pathlib import Path
from typing import Union, Dict, Any
import numpy as np
from pandas.core.frame import DataFrame
import seaborn as sns
import pandas as pd
import copy
import matplotlib.pyplot as plt

from mutagenesis_visualization.main.classes.base_model import Pyplot


class MultipleKernel(Pyplot):
    """
    Class to generate plots of multiple kernels.
    """
    def __call__(
        self,
        dict_entries,
        colors=['k', 'crimson', 'dodgerblue', 'g', 'silver'],
        output_file: Union[None, str, Path] = None,
        **kwargs
    ):
        """
        Generate a kernel density plot for multiple objects passed as a
        dictionary. Can manage either Screen objects or dataframes out
        of the calculate_enrichments function.

        Parameters
        ----------
        dict_entries : dictionary containing dataframes
            Allows for either putting multiple objects as inputs or to
            use dataframes that come out of the calculate_enrichments
            function. If you use an object, you need to say object.dataframe.

        colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
            List of the colors (in order of arguments) that the kernels
            will have.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # hard copy of data
        dict_copy = copy.deepcopy(dict_entries)

        # create figure
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])

        # plot (allows two types of data input)
        for (label, dataset, color) in zip(dict_copy.keys(), dict_copy.values(),
                                           colors[0 : len(dict_copy)]):
            if isinstance(dataset, DataFrame):  # check if input is a dataframe
                # plot objects scores
                self.ax_object = sns.kdeplot(dataset['Score_NaN'], color=color, lw=2, label=label)
            else:  # otherwise assume its an array
                # get rid of stop codons
                dataset.drop('*', errors='ignore', inplace=True)
                dataset = dataset.stack()
                # plot stacked matrix
                self.ax_object = sns.kdeplot(
                    dataset[~np.isnan(dataset)], color=color, lw=2, label=label
                )

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.xlabel(
            r'$∆E^i_x$',
            fontsize=temp_kwargs["x_label_fontsize"],
            fontname='Arial',
            color='k',
            labelpad=0
        )
        plt.ylabel(
            'Probability density',
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
        plt.legend(
            self.dict_copy.keys(),
            loc='best',
            frameon=False,
            fontsize=9,
            handlelength=1,
            handletextpad=0.5
        )

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
        temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))
        return temp_kwargs
