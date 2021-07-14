import numpy as np
import seaborn as sns
import pandas as pd
import copy
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Union

from mutagenesis_visualization.main.classes.base_model import Pyplot


class MultipleKernel(Pyplot):
    """
    Class to generate
    """
    def __call__(
        dict_entries,
        colors=['k', 'crimson', 'dodgerblue', 'g', 'silver'],
        output_file: Union[None, str, Path] = None,
        **kwargs
    ):
        """
        Generate a kernel density plot for multiple objects passed as a dictionary.
        If specified it can also draw a histogram. Uses sns.dispplot. Can manage either
        Screen objects or dataframes out of the calculate_enrichments function.

        Parameters
        ----------
        dict_entries : dictionary containing dataframes
            Allows for either putting multiple objects as inputs or to use dataframes
            that come out of the calculate_enrichments function. If you use an object,
            you need to say object.dataframe.

        colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
            List of the colors (in order of arguments) that the kernels will have.

        output_file : str, default None
            If you want to export the generated graph, add the path and name of the file.
            Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
            return_plot_object : boolean, default False
                If true, will return plotting objects (ie. fig, ax).

        Returns
        ----------
        fig, ax : matplotlib figure and subplots
            Needs to have return_plot_object=True. By default they do
            not get returned.

        """

        # update kwargs
        temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
        temp_kwargs.update(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
        temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))

        # hard copy of data
        dict_copy = copy.deepcopy(dict_entries)

        # create figure
        fig = plt.figure(figsize=temp_kwargs['figsize'])

        # import parameters
        code_kwargs._parameters()

        # plot (allows two types of data input)
        for (label, dataset, color) in zip(dict_copy.keys(), dict_copy.values(),
                                           colors[0 : len(dict_copy)]):
            if isinstance(dataset, pd.core.frame.DataFrame):  # check if input is a dataframe
                # plot objects scores
                ax = sns.kdeplot(dataset['Score_NaN'], color=color, lw=2, label=label)
            else:  # otherwise assume its an array
                # get rid of stop codons
                dataset.drop('*', errors='ignore', inplace=True)
                dataset = dataset.stack()
                # plot stacked matrix
                ax = sns.kdeplot(dataset[~np.isnan(dataset)], color=color, lw=2, label=label)

        # tune graph
        plt.xlabel(r'$∆E^i_x$', fontsize=10, fontname='Arial', color='k', labelpad=0)
        plt.ylabel('Probability density', fontsize=10, fontname='Arial', color='k', labelpad=3)
        plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
        plt.xlim(temp_kwargs['xscale'])
        plt.grid()
        plt.legend(
            dict_copy.keys(),
            loc='best',
            frameon=False,
            fontsize=9,
            handlelength=1,
            handletextpad=0.5
        )

        # save file
        code_utils._save_work(fig, output_file, temp_kwargs)

        # return matplotlib object
        if temp_kwargs['return_plot_object']:
            return fig, ax

        if temp_kwargs['show']:
            plt.show()
