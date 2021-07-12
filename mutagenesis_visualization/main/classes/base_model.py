"""
This module contains the parent class for all the plot classes.
"""
from pathlib import Path
from typing import Union, Dict, Any, Optional
import copy
from matplotlib import rcParams
from matplotlib.figure import Figure
from mutagenesis_visualization.main.utils.kwargs import default_kwargs


class Pyplot:
    """
    Plot abstract class used to visualize mutagenesis data using
    matplotlib.
    """
    def __init__(self, dataset: Optional[Any] = None) -> None:
        """
        Docstring placeholder
        """
        self.dataset: Any = dataset
        self.kwargs: Dict[str, Any] = default_kwargs()
        self.fig: Figure = Figure()
        self.ax_object: Optional[Any] = None

    def _save_work(self, output_file: Union[None, str, Path], temp_kwargs: Dict[str, Any]) -> None:
        '''
        Save file function using pathlib.
        '''
        if output_file:
            self.fig.savefig(
                Path(output_file),
                format=Path(output_file).suffix.strip('.'),
                bbox_inches='tight',
                dpi=temp_kwargs['dpi'],
                transparent=True
            )

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = copy.deepcopy(self.kwargs)
        return temp_kwargs.update(kwargs)

    def _load_parameters(self, ) -> None:
        """
        Default rcParams.
        """
        # normal font
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Arial']

        # math font
        rcParams['mathtext.fontset'] = 'custom'
        rcParams['mathtext.rm'] = 'Arial'
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
        rcParams['mathtext.fontset'] = 'custom'
        rcParams['mathtext.rm'] = 'Arial'
        rcParams['svg.fonttype'] = 'none'

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        pass

    def return_plot_object(self, ):
        """
        Return matplotlib object.
        """
        return self.fig, self.ax_object
