"""
This module contains the class that plots histograms.
"""
from typing import Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from mutagenesis_visualization.main.classes.base_model import Pyplot


class Histogram(Pyplot):
    """
    Class to generate a histogram plot.
    """
    def __call__(
        self,
        population: str = 'All',
        show_parameters: bool = False,
        loc: str = "best",
        replicate: int = -1,
        output_file: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a histogram plot. Can plot single nucleotide variants
        (SNVs) or non-SNVs only.

        Parameters
        ----------
        population : str, default 'All'.
            Other options are 'SNV' and 'nonSNV'.

        show_parameters: bool, default False
            If set to true, will display the mean and the median of
            the data.

        loc: str, default "best"
            Set the location of the meam and median. Check the matplotlib
            plt.legend method to see how the parameter loc works.
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_file : str, default None
            If you want to export the generated graph, add the path and
            name of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            return_plot_object : boolean, default False
                If true, will return plotting objects (ie. fig, ax).
            bins : int or str, default 'auto'.
                Number of bins for the histogram. By default it will
                automatically decide the number of bins.
            color: str, default 'k'
                Change to a different color if desired.
        """
        temp_kwargs = self._update_kwargs(kwargs)
        self.fig = plt.figure(figsize=temp_kwargs['figsize'])
        self.graph_parameters()

        # Select case input data
        data_to_use = self.dataframes.df_notstopcodons[replicate]['Score_NaN']
        if population == 'SNV':
            data_to_use = self.dataframes.df_snv[replicate]['Score_NaN']
        elif population == 'nonSNV':
            data_to_use = self.dataframes.df_nonsnv[replicate]['Score_NaN']

        # plot histogram
        self.ax_object = plt.hist(data_to_use, density=True, bins=temp_kwargs['bins'], color=temp_kwargs['color'])

        # calculate parameters
        if show_parameters:
            plt.plot(0) # so there are two labels
            legend_labels = ('x̄ = '+str(round(data_to_use.mean(),3)),'x̃ = '+str(round(data_to_use.median(),3)))
            plt.legend(legend_labels,frameon=False, framealpha=1, loc=loc,
                borderaxespad=0, ncol=1,handlelength=0, handletextpad=0)

        self._tune_plot(temp_kwargs)
        self._save_work(output_file, temp_kwargs)

        if temp_kwargs['show']:
            plt.show()

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
        return temp_kwargs

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # axes labels and title
        plt.xlabel(
            r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
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

        # axes limits. spacer will be 1 or the
        if temp_kwargs['xscale'] != (None, None):
            plt.xlim(temp_kwargs['xscale'])
            plt.xticks(
                np.arange(
                    temp_kwargs['xscale'][0],
                    temp_kwargs['xscale'][1] + temp_kwargs['tick_spacing'],
                    temp_kwargs['tick_spacing']
                )
            )
        if temp_kwargs['yscale'] != (None, None):
            plt.ylim(temp_kwargs['yscale'])
            plt.grid()
