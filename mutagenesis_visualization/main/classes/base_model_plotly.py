"""
This module contains the parent class for all the plot classes.
"""
from pathlib import Path
from typing import Union, Dict, Any, Optional, List
from copy import deepcopy
from pandas import DataFrame
from numpy import typing as npt
from plotly.graph_objects import Figure
from mutagenesis_visualization.main.utils.kwargs import generate_default_kwargs
from mutagenesis_visualization.main.utils.replicates_screen_input import DataframesHolder


class Plotly:
    """
    Plot abstract class used to visualize mutagenesis data using
    plotly.
    """
    def __init__(
        self,
        dataframes: DataframesHolder,
        aminoacids: Union[str, List[str]] = "",
        datasets: Optional[List[npt.NDArray]] = None,
        sequence: Optional[str] = None,
        start_position: Optional[int] = None,
        end_position: Optional[int] = None,
    ) -> None:
        """
        Docstring placeholder
        """
        if isinstance(aminoacids, str):
            self.aminoacids: List[str] = list(aminoacids)
        elif isinstance(aminoacids, list):
            self.aminoacids = aminoacids

        if datasets:
            self.datasets: List[npt.NDArray] = datasets
        self.dataframes: DataframesHolder = dataframes
        self.sequence: Optional[str] = sequence
        self.start_position: Optional[int] = start_position
        self.end_position: Optional[int] = end_position
        self.kwargs: Any = generate_default_kwargs()
        self.fig: Any = None
        self.df_output: Optional[DataFrame] = None

    def _save_html(self, output_html: Union[None, str, Path], temp_kwargs: Dict[str, Any]) -> None:
        """
        Save figure to html.
        """
        if output_html:
            Path(output_html).parent.mkdir(parents=True, exist_ok=True)
            if '.html' not in str(output_html):
                output_html = str(output_html) + '.html'
            self.fig.write_html(str(Path(output_html)))
        if temp_kwargs['show']:
            self.fig.show(config={'displayModeBar': False})

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = deepcopy(self.kwargs)
        temp_kwargs.update(kwargs)
        temp_kwargs['aminoacids'] = kwargs.get('aminoacids', self.aminoacids)
        if "*" in self.aminoacids and len(self.aminoacids) == 21:
            temp_kwargs['neworder_aminoacids'] = kwargs.get(
                'neworder_aminoacids', list('DEKHRGNQASTPCVYMILFW*')
            )
        elif len(self.aminoacids) == 20:
            temp_kwargs['neworder_aminoacids'] = kwargs.get(
                'neworder_aminoacids', list('DEKHRGNQASTPCVYMILFW')
            )
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
