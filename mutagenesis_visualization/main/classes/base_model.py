"""
This module contains the parent class for all the plot classes.
"""
from pathlib import Path
from typing import Tuple, Union, Dict, Any, Optional, List
import copy
from matplotlib import rcParams
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.utils.kwargs import generate_default_kwargs


class Pyplot:
    """
    Plot abstract class used to visualize mutagenesis data using
    matplotlib.
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
        sequence_raw: Optional[str] = None,
        start_position: Optional[int] = None,
        secondary: Optional[list] = None,
        secondary_dup: Optional[list] = None,
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
        self.sequence_raw: Optional[str] = sequence_raw
        self.start_position: Optional[int] = start_position
        self.secondary: Optional[list] = secondary
        self.secondary_dup: Optional[list] = secondary_dup
        self.kwargs: Dict[str, Any] = generate_default_kwargs()
        self.fig: Figure = Figure()
        self.ax_object: Optional[Any] = None
        self.cb_object: Optional[Any] = None
        self.gs_object: Optional[Any] = None
        self.df_output: Optional[DataFrame] = None
        self.screen_object: Optional[Any] = None
        self.sequence_updated: Optional[str] = None
        self.ax_object2 = None
        self.ax_object3 = None
        self.average_residue = None

    def _save_work(self, output_file: Union[None, str, Path], temp_kwargs: Dict[str, Any]) -> None:
        """
        Save file function using pathlib.
        """
        if output_file:
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

    def _load_parameters(self, ) -> None:
        """
        Default rcParams.
        """
        # normal font
        #rcParams['font.family'] = 'sans-serif'
        #rcParams['font.sans-serif'] = ['Arial']

        # math font
        #rcParams['mathtext.fontset'] = 'custom'
        #rcParams['mathtext.rm'] = 'Arial'
        rcParams['svg.fonttype'] = 'none'

        # add grid
        rcParams['grid.color'] = 'silver'
        rcParams['grid.linestyle'] = '--'
        rcParams['grid.linewidth'] = 1
        rcParams['lines.dashed_pattern'] = [5, 10]
        rcParams['axes.axisbelow'] = True
        # Parameters for all graphs
        rcParams['xtick.labelsize'] = 9
        rcParams['ytick.labelsize'] = 9

    def _font_parameters(self, ) -> None:
        """
        Default math font rcParams.
        """
        # math font
        #rcParams['mathtext.fontset'] = 'custom'
        #rcParams['mathtext.rm'] = 'Arial'
        rcParams['svg.fonttype'] = 'none'

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
