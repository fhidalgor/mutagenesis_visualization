"""
This module contains the kernel class.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
from seaborn import kdeplot
from mutagenesis_visualization.main.classes.base_model import Pyplot


class Kernel(Pyplot):
    """
    Class to generate a kernel density plot.
    """
    def __call__(
        self,
        show_replicates: bool = False,
        wt_counts_only: bool = False,
        show_mean: bool = False,
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot univariate or bivariate distributions using kernel density estimation.

        Parameters
        ----------
        show_replicates: bool, optional default False
            If set to true, will plot the kernel of each replicate.

        wt_counts_only: bool, optional default False
            If set to true, it will plot the kernel distribution of the
            wild type alleles only.

        show_mean: bool, optional default False
            If set to true, it will plot the kernel distribution mean of
            replicates when wt_counts_only is True. Otherwise, it will
            show the mean by default so this parameter won't work if
            show_replicates is set to False.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of the file.
            Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            color: str, default "k"
                Set the color of the mean plot.

            kernel_color_replicates : list of colors, default None
                Add a list of color codes to tune the colors of the plots.

            return_plot_object : boolean, default False
                If true, will return plotting objects (ie. fig, ax_object).
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self.graph_parameters()

        self.fig = plt.figure(figsize=temp_kwargs['figsize'])

        # plot kernel
        if show_replicates:
            assert len(self.dataframes.df_notstopcodons) > 1, "No replicates found."
            if wt_counts_only:
                data_to_plot = self.dataframes.df_wildtype_scores
            else:
                data_to_plot = self.dataframes.df_notstopcodons

            for i, df_replicate in enumerate(data_to_plot[:-1], start=1):
                label = "r{}".format(str(i))
                self.ax_object = kdeplot(
                    df_replicate['Score_NaN'],
                    color=temp_kwargs["kernel_color_replicates"][i-1],
                    lw=2,
                    label=label
                )
            if show_mean:
                self.ax_object = kdeplot(
                    data_to_plot[-1]['Score_NaN'], color=temp_kwargs["color"], lw=2, label="mean"
                )

            plt.legend(loc='best', frameon=False, handlelength=1, handletextpad=0.5)
        else:
            if wt_counts_only:
                data_to_plot = self.dataframes.df_wildtype_scores
            else:
                data_to_plot = self.dataframes.df_notstopcodons
            self.ax_object = kdeplot(
                data_to_plot[replicate]['Score_NaN'],
                color=temp_kwargs["color"],
                lw=2,
            )
            self.ax_object.get_legend().remove()
        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        plt.xlabel(
            temp_kwargs['x_label'],
            fontsize=temp_kwargs["x_label_fontsize"],
            color='k',
            labelpad=0
        )
        plt.ylabel(
            temp_kwargs['y_label'],
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
        plt.grid()

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
        temp_kwargs['x_label'] = kwargs.get('x_label', r'$∆E^i_x$')
        temp_kwargs['y_label'] = kwargs.get('y_label', 'Probability density')
        temp_kwargs['kernel_color_replicates'] = kwargs.get('kernel_color_replicates', [None] * 50)
        return temp_kwargs
