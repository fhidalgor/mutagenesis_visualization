"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
import numpy as np
from mutagenesis_visualization.main.classes.base_model import Pyplot
#from mutagenesis_visualization.main.utils.u import Counter

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
    Class to generate a ROC analysis.
    """
    def plot(self, output_file: Union[None, str, Path] = None, **kwargs: Dict[str, Any],):
        """
        Generates a bar plot of data sorted by secondary elements (alpha
        helices and beta sheets).

        Parameters
        -----------
        self : object from class *Screen*

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # Get data
        df_output: DataFrame = _calculate_secondary(self.dataframe, self.secondary_dup)

        # Color
        df_output['Color'] = df_output.apply(code_utils._color_data, axis=1, args=(temp_kwargs['color_gof'], temp_kwargs['color_lof'],))

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
        self.ax_object.set_xticklabels(labels, fontsize=9, fontname="Arial", color='k', minor=False, rotation=0,)
        self.ax_object.set_ylabel(r'$∆E^i_x$', fontsize=10, fontname="Arial", color='k', labelpad=12, rotation=0,)
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        plt.title(temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=10, pad=5,)

        self._save_work(output_file, temp_kwargs)

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
        temp_kwargs['yscale'] = kwargs.get('yscale', (-2, 1))
        return temp_kwargs
