"""
This module contains the group correlation class.
"""
from typing import List, Union, Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
#import logomaker

from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.pca_utils import calculate_substitution_correlations


class GroupCorrelation(Pyplot):
    """
    This class will conduct a correlation from the enrichment scores.
    """
    def plot(
        self,
        r2: float = 0.5,
        groups: List[str]=['DEHKR', 'QN', 'CASTG', 'ILMV', 'WYF'],
        output_file: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ):
        """
        Determines which amino acids better represent the heatmap. Requires
        logomaker package.

        Parameters
        -----------
        self : object from class *Screen*

        r2 : float
            cutoff of the r**2 correlation value. Only values above that will
            be plot at the sequence logo.

        groups : list, default ['DEHKR','QN','CASTG','ILMV','WYF']
            groups of aa to combine together

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        self._load_parameters()

        # If there is a stop codon, delete it
        if '*' in temp_kwargs['neworder_aminoacids']:
            temp_kwargs['neworder_aminoacids'].remove('*')

        # Get R2 of each combination of amino acid substitutions
        df = calculate_substitution_correlations(self, temp_kwargs['neworder_aminoacids'], groups)

        # Filter according the the R2 correlation value
        filtered = df.loc[df['R2'] > r2]
        logoplot = logomaker.alignment_to_matrix(list(filtered['Combinations']))

        # create Logo object
        self.fig = logomaker.Logo(
            logoplot,
            font_name='Arial',
            color_scheme='chemistry',
            vpad=.1,
            width=.8,
            figsize=((len(logoplot) + 1) / 2.5, 1)
        )

        # style using Logo methods
        self.fig.style_xticks(anchor=0, spacing=1, rotation=0)

        # No yticks and no xticks (but keep labels)
        plt.yticks([], [])
        self.fig.ax.tick_params(axis='both', which='both', length=0)

        # style using Axes methods
        self.fig.ax.set_ylabel('Bits')
        self.fig.ax.set_xlim([-0.5, len(logoplot) - 0.5])

        # for putting title on graph
        plt.title(
            temp_kwargs['title'],
            horizontalalignment='center',
            fontname="Arial",
            fontsize=10,
            pad=10,
        )

        self._save_work(output_file, temp_kwargs)
