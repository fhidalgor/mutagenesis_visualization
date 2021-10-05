"""
This module contains the plot cumulative enrichment class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from matplotlib import ticker

from mutagenesis_visualization.main.classes.base_model import Pyplot


class Cumulative(Pyplot):
    """
    This class will plot a cumulative function on the enrichment scores
    from first to last amino acid.
    """
    def __call__(
        self,
        mode: str = 'all',
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates a cumulative plot of the enrichment scores by position.

        Parameters
        -----------
        mode : str, default 'all'
            Options are 'mean', 'all','SNV' and 'nonSNV'.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # create figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # Get data filtered
        self.df_output = self._filter_by_snv(replicate, mode)
        cumsum = self.df_output.cumsum(skipna=False)['Score']
        plt.plot(self.df_output['Position'], cumsum / list(cumsum)[-1], color='red', lw=2)

        temp_kwargs['y_label'] = 'Cumulative LoF'
        if list(cumsum)[-1] > 0:
            temp_kwargs['y_label'] = 'Cumulative GoF'

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Graph limits
        plt.xlim(
            self.dataframes.df_notstopcodons[-1]['Position'].min(),
            self.dataframes.df_notstopcodons[-1]['Position'].max() + 1
        )
        plt.ylim(0, 1.1)

        self.ax_object.xaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))

        # Axis labels
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            color='k'
        )
        plt.ylabel(temp_kwargs['y_label'], fontsize=temp_kwargs['y_label_fontsize'], color='k', labelpad=5)
        plt.xlabel('Position', fontsize=temp_kwargs["x_label_fontsize"], color='k', labelpad=0)

        # x=y line
        plt.plot([0, self.df_output['Position'].max()], [0, 1],
                 color='silver',
                 lw=2,
                 linestyle='--')

    def _filter_by_snv(self, replicate: int, mode: str) -> DataFrame:
        if mode.lower() == 'all':
            return self.dataframes.df_notstopcodons[replicate]
        if mode.lower() == 'snv':
            return self.dataframes.df_snv[replicate]
        if mode.lower() == 'nonsnv':
            return self.dataframes.df_nonsnv[replicate]
        if mode.lower() == 'mean':
            return self.dataframes.df_notstopcodons[replicate].groupby(
                by='Position', as_index=False
            ).mean()
        raise ValueError("mode not selected")

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        # load labels
        temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2))
        temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 20)
        return temp_kwargs
