"""
This module has the class for a bar plot of the library representation.
"""
from pathlib import Path
from typing import Union, Dict, Any, List, Optional
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker

from mutagenesis_visualization.main.classes.base_model import Pyplot


def _percentage_column(df_input: pd.DataFrame):
    """
    Make the percentage per column.
    """
    return df_input / df_input.sum(axis=0) * 100


class LibraryRepresentation(Pyplot):
    """
    Class to generate a library representation bar plot.
    """
    def __init__(
        self,
        dataset: pd.DataFrame,
        aminoacids: List[str],
        positions: List[int],
    ) -> None:
        super().__init__(dataset=dataset)
        self.aminoacids: List[str] = aminoacids
        self.positions: List[int] = positions
        self.df_percentage: Optional[pd.DataFrame] = None

    def plot(
        self,
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Generates a cumulative stacked bar plot. Each bar represents an
        amino acid position, and each color indicates the observed
        variant frequency.

        Parameters
        ----------
        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
         """
        temp_kwargs = self._update_kwargs(kwargs)
        self._load_parameters()

        # Transform data
        self.df_percentage = _percentage_column(self._group_codons_to_aa(self.dataset))

        # colors
        colors = plt.cm.tab20(np.linspace(0, 1, 20))

        # plot
        self.fig, ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        self.ax_object = self.df_percentage.T.plot.bar(
            stacked=True,
            ax=ax_object,
            color=colors,
        )
        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _tune_plot(self, temp_kwargs) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes parameters
        self.ax_object.set_ylim(0, 100)
        self.ax_object.set_ylabel(
            'Cumulative % AA representation',
            fontsize=12,
            fontname="Arial",
            color='k',
            labelpad=0,
            rotation=90,
        )
        self.ax_object.yaxis.set_major_formatter(ticker.PercentFormatter())

        self.ax_object.set_xlabel(
            'Amino acid position',
            fontsize=12,
            fontname="Arial",
            color='k',
            labelpad=4,
        )
        self.ax_object.set_xticklabels(self.positions)

        plt.title(temp_kwargs['title'], fontsize=14, fontname='Arial', color='k', pad=20)

        # legend
        plt.legend(
            markerscale=0.5,
            frameon=False,
            framealpha=0,
            handlelength=0.75,
            handletextpad=0.25,
            ncol=len(self.df_percentage),
            columnspacing=0.75,
            bbox_to_anchor=(0.5, 1.09),
            loc='upper center'
        )

    def _update_kwargs(self, kwargs) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (10, 4))
        return temp_kwargs

    def _group_codons_to_aa(self, df_input: pd.DataFrame):
        """
        Group different codons that are synonymous. Returns sum of counts.
        """
        df_input = df_input.copy()
        df_input['Aminoacid'] = self.aminoacids
        # Group by mean
        return df_input.groupby(as_index=True, by='Aminoacid', sort=False).sum()
