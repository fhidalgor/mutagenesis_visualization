from typing import Union, Dict, Any
from pathlib import Path
import copy
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from matplotlib import ticker

from mutagenesis_visualization.main.classes.base_model import Pyplot

class Cumulative(Pyplot):
    """
    This class will plot a cumulative function on the enrichment scores
    from first to last amino acid.
    """
    def __init__(
        self,
        dataframe: DataFrame,
        dataframe_snv: DataFrame,
        dataframe_nonsnv: DataFrame,
    ) -> None:
        super().__init__(
            dataframe=dataframe,
        )
        self.dataframe_snv: DataFrame = dataframe_snv
        self.dataframe_nonsnv: DataFrame=dataframe_nonsnv

    def plot(self, mode: str ='all', output_file: Union[None, str, Path] = None, **kwargs: Dict[str, Any],):
        """
        Generates a cumulative plot of the enrichment scores by position.

        Parameters
        -----------
        self : object from class *Screen*

        mode : str, default 'all'
            Options are 'mean', 'all','SNV' and 'nonSNV'.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # create figure
        self.fig, self.ax_object = plt.subplots(figsize=temp_kwargs['figsize'])

        # Get data filtered
        self.df_output = self._filter_by_snv(mode)
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
        plt.xlim(self.dataframe['Position'].min(), self.dataframe['Position'].max() + 1)
        plt.ylim(0, 1.1)

        self.ax_object.xaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))

        # Axis labels
        plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
        plt.ylabel(temp_kwargs['y_label'], fontsize=12, fontname="Arial", color='k', labelpad=5)
        plt.xlabel('Position', fontsize=12, fontname="Arial", color='k', labelpad=0)

        # x=y line
        plt.plot([0, self.df_output['Position'].max()], [0, 1], color='silver', lw=2, linestyle='--')

    def _filter_by_snv(self, mode: str) -> DataFrame:
        if mode.lower() == 'all':
            return self.dataframe
        elif mode.lower() == 'snv':
            return self.dataframe_snv
        elif mode.lower() == 'nonsnv':
            return self.dataframe_nonsnv
        elif mode.lower() == 'mean':
            return self.dataframe.groupby(by='Position', as_index=False).mean()

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        # load labels
        temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2))
        temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 20)
        return temp_kwargs
