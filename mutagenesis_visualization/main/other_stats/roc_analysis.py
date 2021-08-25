"""
This module contains the box plot class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from matplotlib import ticker
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.other_stats_utils import (
    select_grouping,
    merge_class_variants,
    roc_auc,
)


class ROC(Pyplot):
    """
    Class to generate a ROC analysis.
    """
    def __call__(
        self,
        df_class: DataFrame,
        mode: str = 'pointmutant',
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates ROC AUC plot. It compares enrichment scores to some labels
        that the user has specified.

        Parameters
        -----------
        df_class: Pandas dataframe
            A dataframe that contains a column of variants labeled 'Variant'
            with a column labeled 'Class' containing the true class of that
            mutation.

        mode : str, default 'pointmutant'
            Specify what enrichment scores to show. If mode = 'mean', it will
            show the mean of each position. If mode = 'A', it will show the
            alanine substitution profile. Can be used for each amino acid.
            Use the one-letter code and upper case.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # Chose mode: !!! not used
        df_grouped: DataFrame = select_grouping(self.dataframe, mode)

        # Merge dataframe with classes
        df_output: DataFrame = merge_class_variants(self.dataframe, df_class, mode)

        # Calculate ROC parameters
        fpr, tpr, auc, _ = roc_auc(df_output)

        # create figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        lw: int = 2
        plt.plot(fpr, tpr, color='k', lw=lw, label='AUC = %0.2f' % auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Graph limits
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        tick_spacing: float = 0.2
        self.ax_object.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        self.ax_object.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        # Axis labels
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            fontname='Arial',
            color='k'
        )
        plt.ylabel(
            'True Positive Rate',
            fontsize=12,
            fontname="Arial",
            color='k',
            labelpad=0,
        )
        plt.xlabel('False Positive Rate', fontsize=12, fontname="Arial", color='k')

        # Legend
        plt.legend(
            loc='lower right',
            handlelength=0,
            handletextpad=0,
            frameon=False,
        )

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2.5))
        return temp_kwargs