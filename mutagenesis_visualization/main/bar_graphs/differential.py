"""
This module contains the class that plots the differential bar plot.
"""
from typing import Union, Dict, Any
from pathlib import Path
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import process_rmse_residue
from mutagenesis_visualization.main.utils.heatmap_utils import generate_cartoon


class Differential(Pyplot):
    """
    Class to generate the difference between two experiments.
    """
    def __call__(
        self,
        screen_object: Any,
        metric: str = 'rmse',
        plot_type: str = 'bar',
        show_cartoon: bool = False,
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Plot the mean positional difference between two experiments.

        Parameters
        ----------
        screen_object : another *Screen* object to compare with.

        metric: str, default 'rmse'
            The way to compare the two objects.
            Options are 'rmse' ((x-y)**2/N)**0.5, 'squared' ((x**2-y**2)/N and
            'mean' (x-y)/N.

        plot_type: str, default 'bar'
            Options are 'bar' and 'line'.

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
        kwargs["metric"] = metric.lower()
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        self.df_output: DataFrame = process_rmse_residue(
            self.dataframe,
            screen_object.dataframe,
            temp_kwargs['metric']
        )
        # make cartoon
        if show_cartoon:
            self.fig = plt.figure(figsize=temp_kwargs['figsize'])
            self.gs_object = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
            self.ax_object = plt.subplot(self.gs_object[0])
        else:
            self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # plot
        if plot_type.lower() == 'line':
            self.ax_object.plot(self.df_output['Position'], self.df_output['d1 - d2'], color=temp_kwargs['color'])
        else:
            self.ax_object.bar(
                self.df_output['Position'],
                self.df_output['d1 - d2'],
                width=1.2,
                color=temp_kwargs['color'],
                snap=False
            )

        # Needs to be worked for selecting the longer one
        # cartoon
        title_pad = 0
        if show_cartoon:
            title_pad = 2.5
            generate_cartoon(
                self.secondary,
                self.start_position,
                self.gs_object,
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
        temp_kwargs: Dict[str, Any] =  super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2.5))
        temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 5)
        temp_kwargs['color'] = kwargs.get('color', 'grey')
        temp_kwargs['x_label'] = kwargs.get('x_label', 'Position')

        if kwargs["metric"] == 'mean':
            temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
            temp_kwargs['y_label'] = kwargs.get('y_label', r'Mean Differential $∆E^i_x$')
        elif kwargs["metric"] == 'rmse':
            temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
            temp_kwargs['y_label'] = kwargs.get('y_label', r'RMSE Differential')
        if kwargs["metric"] == 'squared':
            temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
            temp_kwargs['y_label'] = kwargs.get('y_label', r'Squared Differential')
        return temp_kwargs

    def _tune_plot(self, temp_kwargs) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # add grid
        plt.grid()

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
                      len(self.df_output) + self.start_position, temp_kwargs['tick_spacing'])
        )
        self.ax_object.set_xlabel('Residue', fontsize=10, fontname="Arial", color='k', labelpad=4)
        self.ax_object.set_xlim(
            self.start_position - 0.1,
            len(self.df_output) + self.start_position - 1 + 0.1
        )
