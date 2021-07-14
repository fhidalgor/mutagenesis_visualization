"""
This module contains demo data that can be loaded.
"""
from typing import Dict, Union
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.utils.pandas_functions import parse_pivot

PDB_5P21 = "mutagenesis_visualization/data/5p21.pdb"
PDB_1ERM = "mutagenesis_visualization/data/1erm.pdb"
PDB_1A5R = "mutagenesis_visualization/data/1a5r.pdb"
PDB_1ND4 = "mutagenesis_visualization/data/1nd4.pdb"
DEMO_FASTA = "mutagenesis_visualization/data/Ras_family_trimmed.fasta"


def demo(figure: str = 'heatmap', show: bool = True) -> None:
    """
    Performs a demonstration of the mutagenesis_visualization software.

    Parameters
    -----------
    figure : str, default 'heatmap'
        There are a few example plots that can be displayed to test the
        package is working on your station. The options are 'heatmap',
        'miniheatmap', 'mean', 'kernel', 'pca' 'position', 'secondary_mean',
        'correlation', 'individual_correlation'.
        Check the documentation for more information.

    show : boolean, default True
        If True, will do plt.show() for each figure.
    """

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt(
        "mutagenesis_visualization/data/HRas166_RBD.csv", delimiter=','
    )

    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

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
    hras_RBD: Screen = Screen(
        dataset=hras_enrichment_RBD, sequence=hras_sequence, secondary=secondary
    )

    if figure == 'heatmap':
        # Create heatmap plot
        hras_RBD.heatmap(title='H-Ras 2-166', show_cartoon=True, show=show)
    elif figure == 'miniheatmap':
        # Condensed heatmap
        hras_RBD.miniheatmap(title='Wt residue H-Ras', show=show)
    elif figure == 'mean':
        # Mean enrichment by position
        hras_RBD.mean(
            figsize=[6, 2.5], mode='mean', show_cartoon=True, yscale=[-2, 0.5], title='', show=show
        )
    elif figure == 'kernel':
        # Plot kernel dist using sns.distplot.
        hras_RBD.kernel(histogram=True, title='H-Ras 2-166', xscale=[-2, 1], show=show)
    elif figure == 'pca':
        # PCA by amino acid substitution
        hras_RBD.pca(dimensions=[0, 1], figsize=(2, 2), adjustlabels=True, title='', show=show)
    elif figure == 'position':
        # Create plot for position 117
        hras_RBD.position(
            position=117,
            yscale=(-1.5, 0.8),
            figsize=(3.5, 2),
            title='Position 117',
            output_file=None,
            show=show,
        )
    elif figure == 'secondary_mean':
        hras_RBD.secondary_mean(
            yscale=[-1, 0],
            figsize=[3, 2],
            title='Mean of secondary motifs',
            show=show,
        )
    elif figure == 'correlation':
        # Correlation between amino acids
        hras_RBD.correlation(colorbar_scale=[0.5, 1], title='Correlation', show=show)
    elif figure == 'individual_correlation':
        # Explained variability by amino acid
        hras_RBD.individual_correlation(
            yscale=[0, 0.6],
            title='Explained variability by amino acid',
            output_file=None,
            show=show,
        )
    else:
        raise NameError('Select a valid name for a demo figure.')
    return


def demo_datasets() -> Dict[str, DataFrame]:
    """
    Loads example datasets so the user can play with it.

    Returns
    --------
    data_dict : Dict[str, DataFrame]
        Dictionary that contains the datasets used to create the plots on the documentation.

    """
    # Create dictionary where to store data
    data_dict: Dict[str, Union[np.array, DataFrame]] = {}

    # Retrieve H-Ras datasets and store in dict
    hras_enrichment_RBD = np.genfromtxt(
        "mutagenesis_visualization/data/HRas166_RBD.csv",
        delimiter=',',
    )
    data_dict['array_hras_RBD'] = hras_enrichment_RBD

    hras_enrichment_GAPGEF = np.genfromtxt(
        "mutagenesis_visualization/data/HRas166_GAPGEF.csv",
        delimiter=',',
    )
    data_dict['array_hras_GAPGEF'] = hras_enrichment_GAPGEF

    # Beta lactamase data
    df_bla_raw = pd.read_pickle("mutagenesis_visualization/data/df_bla_raw.pkl")
    data_dict['df_bla'], sequence_bla = parse_pivot(df_bla_raw, col_data='DMS_amp_625_(b)')

    # Sumo
    df_sumo1_raw = pd.read_pickle("mutagenesis_visualization/data/df_sumo1_raw.pkl")
    data_dict['df_sumo1'], sequence_sumo1 = parse_pivot(df_sumo1_raw, col_data='DMS')

    # MAPK1
    df_mapk1_raw = pd.read_pickle("mutagenesis_visualization/data/df_mapk1_raw.pkl")
    data_dict['df_mapk1'], sequence_mapk1 = parse_pivot(df_mapk1_raw, col_data='DMS_DOX')

    # UBE2I
    df_ube2i_raw = pd.read_pickle("mutagenesis_visualization/data/df_ube2i_raw.pkl")
    data_dict['df_ube2i'], sequence_ube2i = parse_pivot(df_ube2i_raw, col_data='DMS')

    # TAT
    data_dict['df_tat'] = pd.read_pickle("mutagenesis_visualization/data/df_tat.pkl")

    # REV
    data_dict['df_rev'] = pd.read_pickle("mutagenesis_visualization/data/df_rev.pkl")

    # asynuclein
    data_dict['df_asynuclein'] = pd.read_pickle("mutagenesis_visualization/data/df_asynuclein.pkl")

    # APH
    data_dict['df_aph'] = pd.read_pickle("mutagenesis_visualization/data/df_aph.pkl")

    # b11L5
    df_b11L5F_raw = pd.read_pickle("mutagenesis_visualization/data/df_b11L5F_raw.pkl")
    data_dict['df_b11L5F'], sequence_b11L5F = parse_pivot(
        df_b11L5F_raw, col_data='relative_tryp_stability_score'
    )

    return data_dict
