"""
This module contains the parent class for all the plot classes.
"""
from pathlib import Path
from typing import Tuple, Union, Dict, Any, Optional, List, TYPE_CHECKING
from copy import deepcopy
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from numpy import typing as npt
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.utils.kwargs import generate_default_kwargs
from mutagenesis_visualization.main.utils.matplotlib_parameters import (
    graph_parameters,
    font_parameters,
)
from mutagenesis_visualization.main.utils.replicates_screen_input import DataframesHolder

if TYPE_CHECKING:
    from mutagenesis_visualization.main.classes.screen import Screen


class Pyplot:
    """
    Plot abstract class used to visualize mutagenesis data using
    matplotlib.
    """
    def __init__(
        self,
        dataframes: Optional[DataframesHolder] = None,
        dataframes_raw: Optional[List[DataFrame]] = None,
        aminoacids: Union[str, List[str]] = "",
        datasets: Optional[List[npt.NDArray]] = None,
        sequence: str = "",
        sequence_raw: str = "",
        start_position: int = 0,
        secondary: Optional[List] = None,
        secondary_dup: Optional[List] = None,
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
        if dataframes:
            self.dataframes: DataframesHolder = dataframes
        if dataframes_raw:
            self.dataframes_raw: List[DataFrame] = dataframes_raw
        self.sequence: str = sequence
        self.sequence_raw: str = sequence_raw
        self.start_position: int = start_position
        if secondary:
            self.secondary: list = secondary
        if secondary_dup:
            self.secondary_dup: list = secondary_dup
        self.kwargs: Any = generate_default_kwargs()
        self.fig: Figure = Figure()
        self.ax_object: plt.Axes = None
        self.cb_object: plt.Axes = None
        self.gs_object: plt.Axes = None
        self.df_output: DataFrame = DataFrame()
        self.screen_object: 'Screen' = None
        self.sequence_updated: Union[str, List[str]] = []
        self.ax_object2: plt.Axes = None
        self.ax_object3: plt.Axes = None
        self.average_residue: plt.Axes = None
        self.graph_parameters = graph_parameters
        self.font_parameters = font_parameters
        self.segment: Tuple[int, int] = None

    def _save_work(self, output_file: Union[None, str, Path], temp_kwargs: Dict[str, Any]) -> None:
        """
        Save file function using pathlib.
        """
        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)

            self.fig.savefig(
                Path(output_file),
                format=Path(output_file).suffix.strip('.'),
                bbox_inches='tight',
                dpi=temp_kwargs['dpi'],
                transparent=True
            )
        if temp_kwargs['show']:
            plt.show()

        if temp_kwargs['close']:
            plt.close()

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

    def return_plot_object(self) -> Tuple:
        """
        Return matplotlib object.
        """
        return self.fig, self.ax_object

    def return_dataframe(self) -> DataFrame:
        """
        Return the dataframe that contains the massaged data.
        """
        return self.df_output

    def export_dataframe_to_csv(self, output_file: Union[str, Path]) -> None:
        """
        Return the dataframe that contains the massaged data.
        """
        self.df_output.to_csv(Path(output_file), index=False)
