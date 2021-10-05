"""
This module has the class for a bar plot of the library representation.
"""
from pathlib import Path
from typing import Union, Dict, Any, List, Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.base_model import Pyplot


def _percentage_column(df_input: DataFrame) -> DataFrame:
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
        dataframes_raw: List[DataFrame],
        positions: List[int],
        aminoacids: List[str],
    ) -> None:  # pylint: disable=super-init-not-called
        super().__init__(dataframes_raw=dataframes_raw, aminoacids=aminoacids)
        self.positions: List[int] = positions
        self.df_percentage: Optional[DataFrame] = None

    def __call__(
        self,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates a cumulative stacked bar plot. Each bar represents an
        amino acid position, and each color indicates the observed
        variant frequency.

        Parameters
        ----------
        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
         """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Transform data
        self.df_percentage = _percentage_column(
            self._group_codons_to_aa(self.dataframes_raw[replicate])
        )

        # colors
        colors = plt.cm.get_cmap('tab20')(np.linspace(0, 1, 20))

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

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes parameters
        self.ax_object.set_ylim(0, 100)
        self.ax_object.set_ylabel(
            'Cumulative % AA representation',
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=0,
            rotation=90,
        )
        self.ax_object.yaxis.set_major_formatter(ticker.PercentFormatter())

        self.ax_object.set_xlabel(
            'Amino acid position',
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            labelpad=4,
        )
        self.ax_object.set_xticklabels(self.positions)

        plt.title(temp_kwargs['title'], fontsize=temp_kwargs['title_fontsize'], fontname='Arial', color='k', pad=20)

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

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (10, 4))
        temp_kwargs['title_fontsize'] = kwargs.get('title_fontsize', 14)
        temp_kwargs["y_label_fontsize"] = kwargs.get('y_label_fontsize', 12)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 12)
        return temp_kwargs

    def _group_codons_to_aa(self, df_input: DataFrame) -> DataFrame:
        """
        Group different codons that are synonymous. Returns sum of counts.
        """
        df_input = df_input.copy()
        df_input['Aminoacid'] = self.aminoacids
        # Group by mean
        return df_input.groupby(as_index=True, by='Aminoacid', sort=False).sum()
