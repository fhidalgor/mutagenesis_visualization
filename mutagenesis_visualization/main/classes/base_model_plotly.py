"""
This module contains the parent class for all the plot classes.
"""
from pathlib import Path
from typing import Union, Dict, Any, Optional, List
import copy
from pandas import DataFrame
from plotly.graph_objects import Figure
from mutagenesis_visualization.main.utils.kwargs import generate_default_kwargs


class Plotly:
    """
    Plot abstract class used to visualize mutagenesis data using
    plotly.
    """
    def __init__(
        self,
        aminoacids: Union[str, List[str], None] = None,
        dataset: Optional[Any] = None,
        dataframe: Optional[DataFrame] = None,
        dataframe_stopcodons: Optional[DataFrame] = None,
        dataframe_snv: Optional[DataFrame] = None,
        dataframe_nonsnv: Optional[DataFrame] = None,
        sequence: Optional[str] = None,
        start_position: Optional[int] = None,
        end_position: Optional[int] = None,
    ) -> None:
        """
        Docstring placeholder
        """
        if isinstance(aminoacids, str):
            self.aminoacids: List[str]  = list(aminoacids)
        elif isinstance(aminoacids, list):
            self.aminoacids = aminoacids

        self.dataset: Any = dataset
        self.dataframe: Optional[DataFrame] = dataframe
        self.dataframe_stopcodons: Optional[DataFrame] = dataframe_stopcodons
        self.dataframe_snv: Optional[DataFrame] = dataframe_snv
        self.dataframe_nonsnv: Optional[DataFrame] = dataframe_nonsnv
        self.sequence: Optional[str] = sequence
        self.start_position: Optional[int] = start_position
        self.end_position: Optional[int] = end_position
        self.kwargs: Dict[str, Any] = generate_default_kwargs()
        self.fig: Any = None
        self.df_output: Optional[DataFrame] = None

    def _save_html(self, output_html: Union[None, str, Path], temp_kwargs: Dict[str, Any]) -> None:
        """
        Save figure to html.
        """
        if output_html:
            if '.html' not in output_html:
                output_html += '.html'
            self.fig.write_html(str(Path(output_html)))
        if temp_kwargs['show']:
            self.fig.show(config={'displayModeBar': False})

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['aminoacids'] = kwargs.get('aminoacids', self.aminoacids)
        if "*" in self.aminoacids and len(self.aminoacids)==21:
            temp_kwargs['neworder_aminoacids'] = kwargs.get('neworder_aminoacids', list('DEKHRGNQASTPCVYMILFW*'))
        elif len(self.aminoacids)==20:
            temp_kwargs['neworder_aminoacids'] = kwargs.get('neworder_aminoacids', list('DEKHRGNQASTPCVYMILFW'))
        else:
            temp_kwargs['neworder_aminoacids'] = kwargs.get('neworder_aminoacids', self.aminoacids)
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        pass

    def return_plot_object(self) -> Figure:
        """
        Return plotly object.
        """
        return self.fig
