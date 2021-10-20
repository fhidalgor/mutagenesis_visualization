"""
This module contains the multiple kernel class.
"""

from pathlib import Path
from typing import List, Union, Dict, Any, TYPE_CHECKING, Optional
from seaborn import kdeplot
import matplotlib.pyplot as plt

from mutagenesis_visualization.main.classes.base_model import Pyplot

if TYPE_CHECKING:
    from mutagenesis_visualization.main.classes.screen import Screen


class MultipleKernel(Pyplot):
    """
    Class to generate plots of multiple kernels.
    """
    def __call__(
        self,
        screen_object: Union['Screen', List['Screen']],
        label_kernels: List[str],
        kernel_colors: Optional[List[str]] = None,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any
    ) -> None:
        """
        Generate a kernel density plot for multiple objects.

        Parameters
        ----------
        screen_object : *Screen* object or list containing *Screen*
            objects.

        kernel_colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
            List of the colors (in order of arguments) that the kernels
            will have.

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
        temp_kwargs['label_kernels'] = label_kernels
        self.graph_parameters()

        if not kernel_colors:
            kernel_colors = ['k', 'crimson', 'dodgerblue', 'g', 'silver']

        # create figure
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])

        self.ax_object = kdeplot(
            self.dataframes.df_notstopcodons[replicate]['Score_NaN'],
            color=kernel_colors[0],
            lw=2,
            label=label_kernels[0]
        )
        if isinstance(screen_object, list):
            for (label, sobj, color) in zip(label_kernels, screen_object, kernel_colors[1:]):
                self.ax_object = kdeplot(
                    sobj.dataframes.df_notstopcodons[replicate]['Score_NaN'],
                    color=color,
                    lw=2,
                    label=label
                )
        else:
            self.ax_object = kdeplot(
                screen_object.dataframes.df_notstopcodons[replicate]['Score_NaN'],
                color=kernel_colors[1],
                lw=2,
                label=label_kernels[1]
            )

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.xlabel(
            r'$∆E^i_x$',
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            labelpad=0
        )
        plt.ylabel(
            'Probability density',
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=3
        )
        plt.title(
            temp_kwargs['title'],
            fontsize=temp_kwargs["title_fontsize"],
            color='k'
        )
        plt.xlim(temp_kwargs['xscale'])
        plt.ylim(temp_kwargs['yscale'])
        plt.grid()
        plt.legend(
            temp_kwargs['label_kernels'],
            loc='best',
            frameon=False,
            fontsize=9,
            handlelength=1,
            handletextpad=0.5
        )

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
        return temp_kwargs
