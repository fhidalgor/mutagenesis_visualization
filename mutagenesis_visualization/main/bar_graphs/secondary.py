"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
import numpy as np
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pandas_functions import color_data


def _calculate_secondary(df_input: DataFrame, secondary: list) -> DataFrame:
    """
    Returns copy.
    """
    df_output: DataFrame = df_input.copy()
    df_output.insert(4, 'Secondary', secondary)
    df_output = df_output.groupby(['Secondary'], as_index=False, sort=False).mean()
    df_output = df_output[df_output['Secondary'].str.startswith(('β', 'α'))]
    return df_output.drop(['Position'], axis=1)


class Secondary(Pyplot):
    """
    Class to generate bar plot of data sorted by secondary elements.
    """
    def __call__(
        self,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates a bar plot of data sorted by secondary elements (alpha
        helices and beta sheets).

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
        df_output: DataFrame = _calculate_secondary(
            self.dataframes.df_notstopcodons[replicate], self.secondary_dup
        )

        # Color
        df_output['Color'] = df_output.apply(
            color_data, axis=1, args=(temp_kwargs['color_gof'], temp_kwargs['color_lof'])
        )

        # Make figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        ticks = np.arange(0, len(df_output))  # label locations
        width = 0.5
        labels: list = df_output['Secondary']

        # Plot figure
        self.ax_object.bar(
            ticks,
            df_output['Score'],
            width,
            color=df_output['Color'],
            ec='k',
        )

        # graph parameters
        self.ax_object.set_xticks(ticks)
        self.ax_object.set_xticklabels(
            labels,
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            minor=False,
            rotation=0,
        )
        self.ax_object.set_ylabel(
            r'$∆E^i_x$',
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=12,
            rotation=0,
        )
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontsize=temp_kwargs['title_fontsize'],
            pad=5,
        )

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4, 2))
        temp_kwargs['y_label_fontsize'] = kwargs.get('y_label_fontsize', 10)
        temp_kwargs['title_fontsize'] = kwargs.get('title_fontsize', 10)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 9)
        return temp_kwargs
