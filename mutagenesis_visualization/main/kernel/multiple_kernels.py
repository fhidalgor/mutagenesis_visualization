"""
This module contains the multiple kernel class.
"""

from pathlib import Path
from typing import List, Union, Dict, Any
import seaborn as sns
import copy
import matplotlib.pyplot as plt

from mutagenesis_visualization.main.classes.base_model import Pyplot


class MultipleKernel(Pyplot):
    """
    Class to generate plots of multiple kernels.
    """
    def __call__(
        self,
        screen_object: Union[Any,list],
        label_kernels: List[str],
        colors: List[str]=['k', 'crimson', 'dodgerblue', 'g', 'silver'],
        output_file: Union[None, str, Path] = None,
        **kwargs
    ):
        """
        Generate a kernel density plot for multiple objects.

        Parameters
        ----------
        screen_object : *Screen* object or list containing *Screen*
            objects.

        colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
            List of the colors (in order of arguments) that the kernels
            will have.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or
            'path/filename.svg'.

        **kwargs : other keyword arguments
        """

        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        temp_kwargs['label_kernels'] = label_kernels
        self._load_parameters()

        # create figure
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])

        self.ax_object = sns.kdeplot(self.dataframe['Score_NaN'], color=colors[0], lw=2, label=label_kernels[0])
        if isinstance(screen_object, list):
            for (label, sobj, color) in zip(label_kernels, screen_object, colors):
                self.ax_object = sns.kdeplot(sobj.dataframe['Score_NaN'], color=color, lw=2, label=label)
        else:
            self.ax_object = sns.kdeplot(screen_object.dataframe['Score_NaN'], color=colors[1], lw=2, label=label_kernels[1])

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.xlabel(r'$∆E^i_x$', fontsize=10, fontname='Arial', color='k', labelpad=0)
        plt.ylabel('Probability density', fontsize=10, fontname='Arial', color='k', labelpad=3)
        plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
        plt.xlim(temp_kwargs['xscale'])
        plt.grid()
        plt.legend(
            temp_kwargs['label_kernels'],
            loc='best',
            frameon=False,
            fontsize=9,
            handlelength=1,
            handletextpad=0.5
        )
    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] =  super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
        temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))
        return temp_kwargs
