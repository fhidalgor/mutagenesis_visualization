"""
This module has the class for a bar plot of the mean counts.
"""
from pathlib import Path
from typing import Union, Dict, Any, List
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.base_model import Pyplot


def input_text(text_entries: List[str]) -> None:
    """
    The user can input text as a variable by manually giving the coordinates.
    """
    if text_entries:
        for entry in text_entries:
            plt.text(entry[0], entry[1], entry[2])


class MeanCounts(Pyplot):
    """
    Class to generate a mean counts bar plot.
    """
    def __init__(
        self,
        dataframes_raw: List[DataFrame],
        positions: List[int],
        aminoacids: List[str],
    ) -> None:
        super().__init__(dataframes_raw=dataframes_raw, aminoacids=aminoacids)
        self.positions: List[int] = positions

    def __call__(
        self,
        normalize: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot in a bargraph the mean counts for each residue of the protein.

        Parameters
        ----------
        normalize : bool, default False
            If set to true, the mean counts will be normalized so the
            highest value is 100.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            text_labels : list of lists, default empty
                If you want to add a label to the graph, add the coordinates
                and the text. Example: text_labels = [[x0,y0,text0]
                [x1,y1,text1]].
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        # plot
        self.fig, ax_object = plt.subplots(figsize=temp_kwargs['figsize'])
        self.df_output = self.dataframes_raw[replicate].mean()
        if normalize:
            self.df_output = self.df_output*100/self.df_output.max()
        self.ax_object = self.df_output.T.plot.bar(ax=ax_object, color=temp_kwargs['color'])

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes parameters
        self.ax_object.set_ylim(temp_kwargs['yscale'])
        self.ax_object.set_ylabel(
            temp_kwargs['y_label'],
            fontsize=temp_kwargs["y_label_fontsize"],
            color='k',
            labelpad=0,
            rotation=90,
        )
        self.ax_object.set_xlabel(
            'Amino acid position',
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            labelpad=4,
        )
        self.ax_object.set_xticklabels(self.positions)

        # Title
        plt.title(temp_kwargs['title'], fontsize=temp_kwargs['title_fontsize'], color='k', pad=0)

        # Put text labels
        input_text(temp_kwargs['text_labels'])

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (10, 4))
        temp_kwargs['y_label'] = kwargs.get('y_label', 'Mean counts')
        temp_kwargs['title_fontsize'] = kwargs.get('title_fontsize', 14)
        temp_kwargs["y_label_fontsize"] = kwargs.get('y_label_fontsize', 12)
        temp_kwargs["x_label_fontsize"] = kwargs.get('x_label_fontsize', 12)
        temp_kwargs['color'] = kwargs.get('color', "red")

        return temp_kwargs
