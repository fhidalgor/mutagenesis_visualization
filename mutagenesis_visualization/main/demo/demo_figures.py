"""
This module contains demo data that can be loaded.
"""
from typing import List, Literal
import numpy as np
from numpy import typing as npt

from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.utils.data_paths import HRAS_RBD_COUNTS_CSV, PDB_5P21

FIGURE_OPTIONS = Literal['heatmap', 'miniheatmap', 'mean', 'kernel', 'pca',  # pylint: disable=invalid-name
                         'position', 'secondary_mean', 'correlation', 'individual_correlation',
                         'pymol']


def run_demo(figure: FIGURE_OPTIONS = 'heatmap', show: bool = True) -> None:
    """
    Performs a demonstration of the mutagenesis_visualization software.

    Parameters
    -----------
    figure : str, default 'heatmap'
        There are a few example plots that can be displayed to test the
        package is working on your station. The options are 'heatmap',
        'miniheatmap', 'mean', 'kernel', 'pca' 'position', 'secondary_mean',
        'correlation', 'individual_correlation' and 'pymol'.
        Check the documentation for more information.

    show : boolean, default True
        If True, will execute plt.show() for each figure.
    """

    # Load enrichment scores
    hras_enrichment_rbd: npt.NDArray = np.genfromtxt(HRAS_RBD_COUNTS_CSV, delimiter=',')

    # Define protein sequence
    hras_sequence: str = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'  # pylint: disable=line-too-long

    # Set aminoacids
    aminoacids: List[str] = list('ACDEFGHIKLMNPQRSTVWY*')

    # Define secondary structure
    secondary: list = [['L0'],
                       ['β1'] * (9 - 1), ['L1'] * (15 - 9), ['α1'] * (25 - 15), ['L2'] * (36 - 25),
                       ['β2'] * (46 - 36), ['L3'] * (48 - 46), ['β3'] * (58 - 48),
                       ['L4'] * (64 - 58), ['α2'] * (74 - 64), ['L5'] * (76 - 74),
                       ['β4'] * (83 - 76), ['L6'] * (86 - 83), ['α3'] * (103 - 86),
                       ['L7'] * (110 - 103), ['β5'] * (116 - 110), ['L8'] * (126 - 116),
                       ['α4'] * (137 - 126), ['L9'] * (140 - 137), ['β6'] * (143 - 140),
                       ['L10'] * (151 - 143), ['α5'] * (172 - 151), ['L11'] * (190 - 172)]

    # Create object
    hras_rbd: Screen = Screen(
        datasets=hras_enrichment_rbd,
        sequence=hras_sequence,
        aminoacids=aminoacids,
        secondary=secondary
    )

    if figure.lower() == 'heatmap':
        # Create heatmap plot
        hras_rbd.heatmap(
            title='H-Ras 2-166', mask_selfsubstitutions=False, show_cartoon=True, show=show
        )
    elif figure.lower() == 'miniheatmap':
        # Condensed heatmap
        hras_rbd.miniheatmap(mask_selfsubstitutions=False, title='Wt residue H-Ras', show=show)
    elif figure.lower() == 'mean':
        # Mean enrichment by position
        hras_rbd.enrichment_bar(
            figsize=[6, 2.5], mode='mean', show_cartoon=True, yscale=[-2, 0.5], title='', show=show
        )
    elif figure.lower() == 'kernel':
        # Plot kernel dist using sns.distplot.
        hras_rbd.kernel(histogram=True, title='H-Ras 2-166', xscale=[-2, 1], show=show)
    elif figure.lower() == 'pca':
        # PCA by amino acid substitution
        hras_rbd.pca(dimensions=(0, 1), figsize=(2, 2), adjustlabels=True, title='', show=show)
    elif figure.lower() == 'position':
        # Create plot for position 117
        hras_rbd.position_bar(
            position=117,
            yscale=(-1.5, 0.8),
            figsize=(3.5, 2),
            title='Position 117',
            output_file=None,
            show=show,
        )
    elif figure.lower() == 'secondary_mean':
        hras_rbd.secondary_mean(
            yscale=[-1, 0],
            figsize=[3, 2],
            title='Mean of secondary motifs',
            show=show,
            show_error_bars=False,
        )
    elif figure.lower() == 'correlation':
        # Correlation between amino acids
        hras_rbd.correlation(colorbar_scale=[0.5, 1], title='Correlation', show=show)
    elif figure.lower() == 'individual_correlation':
        # Explained variability by amino acid
        hras_rbd.individual_correlation(
            yscale=[0, 0.6],
            title='Explained variability by amino acid',
            output_file=None,
            show=show,
        )
    elif figure.lower() == 'pymol':
        hras_rbd.pymol(PDB_5P21)
    else:
        raise NameError('Select a valid name for a demo figure.')
