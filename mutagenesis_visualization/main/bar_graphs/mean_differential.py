"""
This module contains the class that plots the mean differential bar plot.
"""
from typing import Union, Dict, Any
import copy
from pathlib import Path
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import process_mean_residue
from mutagenesis_visualization.main.utils.heatmap_utils import generate_cartoon


class MeanDifferential(Pyplot):
    """
    Class to generate a mean enrichment bar plot.
    """
    def plot(
        self,
        screen_object: Any,
        show_cartoon: bool = False,
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ):
        """
        Plot the mean positional difference between two experiments.

        Parameters
        ----------
        screen_object : another Screen object to compare with

        show_cartoon : boolean, default False
            If true, the plot will display a cartoon with the secondary
            structure. The user must have added the secondary structure
            to the object.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # make pandas
        self.df_output: DataFrame = process_mean_residue(
            self.dataframe,
            self.screen_object.dataframe,
        )

        # make cartoon
        if show_cartoon:
            self.fig = plt.figure(figsize=temp_kwargs['figsize'])
            gs_object = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
            self.ax_object = plt.subplot(gs_object[0])
        else:
            self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # plot
        self.ax_object.plot(self.df_output['Position'], self.df_output['d1 - d2'], color='k')

        # Needs to be worked for selecting the longer one
        # cartoon
        title_pad = 0
        if show_cartoon:
            title_pad = 2.5
            generate_cartoon(
                self.secondary,
                self.start_position,
                gs_object,
                1,
                temp_kwargs['cartoon_colors'],
                bottom_space=-0.78,
                show_labels=False,
            )

        # title
        self.ax_object.set_title(
            temp_kwargs['title'],
            fontsize=12,
            fontname='Arial',
            color='k',
            pad=title_pad,
        )

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _update_kwargs(self, kwargs) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2.5))
        temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
        temp_kwargs['y_label'] = kwargs.get('y_label', r'Mean Differential $∆E^i_x$')
        return temp_kwargs

    def _tune_plot(self, temp_kwargs) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # self.ax_objectes parameters
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        self.ax_object.set_ylabel(
            temp_kwargs['y_label'],
            fontsize=10,
            fontname="Arial",
            color='k',
            labelpad=-5,
            rotation=90
        )
        self.ax_object.set_xticks(
            np.arange(self.start_position,
                      len(self.df_output) + self.start_position, 20)
        )
        self.ax_object.set_xlabel('Residue', fontsize=10, fontname="Arial", color='k', labelpad=4)
        self.ax_object.set_xlim(
            self.start_position - 0.1,
            len(self.df_output) + self.start_position - 1 + 0.1
        )
