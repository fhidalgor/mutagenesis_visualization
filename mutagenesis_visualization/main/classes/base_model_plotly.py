"""
This module contains the parent class for all the plot classes.
"""
from pathlib import Path
from typing import Union, Dict, Any, Optional
import copy
from pandas import DataFrame

from mutagenesis_visualization.main.utils.kwargs import default_kwargs


class Plotly:
    """
    Plot abstract class used to visualize mutagenesis data using
    plotly.
    """
    def __init__(self, dataframe: DataFrame) -> None:
        """
        Docstring placeholder
        """
        self.kwargs: Dict[str, Any] = default_kwargs()
        self.fig: Any = None
        self.dataframe: DataFrame = dataframe
        self.df_output: Optional[DataFrame] = None

    def _save_html(self, output_html: Union[None, str, Path]) -> None:
        '''
        Save figure to html.
        '''
        if output_html:
            self.fig.write_html(str(Path(output_html)))

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        return temp_kwargs.update(kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        pass

    def return_plot_object(self):
        """
        Return matplotlib object.
        """
        return self.fig
