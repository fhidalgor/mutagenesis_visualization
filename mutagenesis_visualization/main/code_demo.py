#!/usr/bin/env python
# coding: utf-8

# # Import modules

# In[ ]:


# Regular libraries
import numpy as np
import pandas as pd
import itertools
from os import path
import os

# Local imports
try:
    from mutagenesis_visualization.main.scripts.code_class import Screen
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    from code_class import Screen
    import code_utils
    __file__ = '__main__'


# # Demo

# In[ ]:


def demo(figure='heatmap', show=True):
    """
    Performs a demonstration of the mutagenesis_visualization software.

    Parameters
    -----------
    figure : str, default 'heatmap'
        There are a few example plots that can be displayed to test the package is working on your station.
        The options are 'heatmap', 'miniheatmap', 'mean', 'kernel', 'pca'
        'position', 'secondary_mean', 'correlation', 'individual_correlation'. 
        Check the documentation for more information.
    
    show : boolean, default True
        If True, will do plt.show() for each figure.
        
    Returns
    -------
    None.
    """
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'HRas166_RBD.csv')
    except NameError:
        my_file = os.path.join('../../data', 'HRas166_RBD.csv')

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')

    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

    # Define secondary structure
    secondary = [['L0'], ['β1'] * (9 - 1), ['L1'] * (15 - 9),
                 ['α1'] * (25 - 15), ['L2'] * (36 - 25), ['β2'] * (46 - 36),
                 ['L3'] * (48 - 46), ['β3'] * (58 - 48), ['L4'] * (64 - 58),
                 ['α2'] * (74 - 64), ['L5'] * (76 - 74), ['β4'] * (83 - 76),
                 ['L6'] * (86 - 83), ['α3'] * (103 - 86), ['L7'] * (110 - 103),
                 ['β5'] * (116 - 110), ['L8'] * (126 - 116),
                 ['α4'] * (137 - 126), ['L9'] * (140 - 137),
                 ['β6'] * (143 - 140), ['L10'] * (151 - 143),
                 ['α5'] * (172 - 151), ['L11'] * (190 - 172)]

    # Create object
    hras_RBD = Screen(
        dataset=hras_enrichment_RBD,
        sequence=hras_sequence,
        secondary=secondary
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
            figsize=[6, 2.5],
            mode='mean',
            show_cartoon=True,
            yscale=[-2, 0.5],
            title='',
            show=show
        )
    elif figure == 'kernel':
        # Plot kernel dist using sns.distplot.
        hras_RBD.kernel(
            histogram=True, title='H-Ras 2-166', xscale=[-2, 1], show=show
        )
    elif figure == 'pca':
        # PCA by amino acid substitution
        hras_RBD.pca(
            dimensions=[0, 1],
            figsize=(2, 2),
            adjustlabels=True,
            title='',
            show=show
        )
    elif figure == 'position':
        # Create plot for position 117
        hras_RBD.position(
            position=117,
            yscale=(-1.5, 0.8),
            figsize=(3.5, 2),
            title='Position 117',
            output_file=None,
            show=show
        )
    elif figure == 'secondary_mean':
        hras_RBD.secondary_mean(
            yscale=[-1, 0],
            figsize=[3, 2],
            title='Mean of secondary motifs',
            show=show
        )
    elif figure == 'correlation':
        # Correlation between amino acids
        hras_RBD.correlation(
            colorbar_scale=[0.5, 1], title='Correlation', show=show
        )
    elif figure == 'individual_correlation':
        # Explained variability by amino acid
        hras_RBD.individual_correlation(
            yscale=[0, 0.6],
            title='Explained variability by amino acid',
            output_file=None,
            show=show
        )
    else:
        raise NameError('Select a valid name for a demo figure')
    return


def demo_datasets():
    """
    Loads example datasets so the user can play with it.

    Parameters
    -----------
    None

    Returns
    --------
    data_dict : dictionary
        Dictionary that contains the datasets used to create the plots on the documentation.

    """

    # Use relative file import to access the data folder
    location = os.path.dirname(os.path.realpath(__file__))

    # Create dictionary where to store data
    data_dict = {}

    # Retrieve H-Ras datasets and store in dict
    my_file = os.path.join(location, '../../data', 'HRas166_RBD.csv')
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')
    data_dict['array_hras_RBD'] = hras_enrichment_RBD

    my_file = os.path.join(location, '../../data', 'HRas166_GAPGEF.csv')
    hras_enrichment_GAPGEF = np.genfromtxt(my_file, delimiter=',')
    data_dict['array_hras_GAPGEF'] = hras_enrichment_GAPGEF

    # Beta lactamase data
    my_file = os.path.join(location, '../../data', 'df_bla_raw.pkl')
    df_bla_raw = pd.read_pickle(my_file)
    data_dict['df_bla'], sequence_bla = code_utils.parse_pivot(
        df_bla_raw, col_data='DMS_amp_625_(b)'
    )

    # Sumo
    my_file = os.path.join(location, '../../data', 'df_sumo1_raw.pkl')
    df_sumo1_raw = pd.read_pickle(my_file)
    data_dict['df_sumo1'], sequence_sumo1 = code_utils.parse_pivot(
        df_sumo1_raw, col_data='DMS'
    )

    # MAPK1
    my_file = os.path.join(location, '../../data', 'df_mapk1_raw.pkl')
    df_mapk1_raw = pd.read_pickle(my_file)
    data_dict['df_mapk1'], sequence_mapk1 = code_utils.parse_pivot(
        df_mapk1_raw, col_data='DMS_DOX'
    )

    # UBE2I
    my_file = os.path.join(location, '../../data', 'df_ube2i_raw.pkl')
    df_ube2i_raw = pd.read_pickle(my_file)
    data_dict['df_ube2i'], sequence_ube2i = code_utils.parse_pivot(
        df_ube2i_raw, col_data='DMS'
    )

    # TAT
    my_file = os.path.join(location, '../../data', 'df_tat.pkl')
    data_dict['df_tat'] = pd.read_pickle(my_file)

    # REV
    my_file = os.path.join(location, '../../data', 'df_rev.pkl')
    data_dict['df_rev'] = pd.read_pickle(my_file)

    # asynuclein
    my_file = os.path.join(location, '../../data', 'df_asynuclein.pkl')
    data_dict['df_asynuclein'] = pd.read_pickle(my_file)

    # APH
    my_file = os.path.join(location, '../../data', 'df_aph.pkl')
    data_dict['df_aph'] = pd.read_pickle(my_file)

    # b11L5
    my_file = os.path.join(location, '../../data', 'df_b11L5F_raw.pkl')
    df_b11L5F_raw = pd.read_pickle(my_file)
    data_dict['df_b11L5F'], sequence_b11L5F = code_utils.parse_pivot(
        df_b11L5F_raw, col_data='relative_tryp_stability_score'
    )

    return data_dict


def demo_pdbs():
    """
    Loads example pdbs so the user can play with it (5p21, 1erm, 1a5r, 1nd4).

    Parameters
    -----------
    None

    Returns
    --------
    data_dict : dictionary
        Dictionary that contains the pdbs used to create the plots on the documentation.

    """

    # Use relative file import to access the data folder
    location = os.path.dirname(os.path.realpath(__file__))

    # Create dictionary where to store data
    pdb_dict = {}

    # H-Ras
    pdb_dict['5p21'] = os.path.join(location, '../../data', '5p21.pdb')

    # Beta lactamase data
    pdb_dict['1erm'] = os.path.join(location, '../../data', '1erm.pdb')

    # Sumo
    pdb_dict['1a5r'] = os.path.join(location, '../../data', '1a5r.pdb')

    # APH
    pdb_dict['1nd4'] = os.path.join(location, '../../data', '1nd4.pdb')

    return pdb_dict


def demo_fasta():
    """
    Loads example fasta so the user can play with it.

    Parameters
    -----------
    None

    Returns
    --------
    data_dict : dictionary
        Dictionary that contains the fasta used to extract the sequence conservation.

    """

    # Use relative file import to access the data folder
    location = os.path.dirname(os.path.realpath(__file__))
    print(location)
    # Create dictionary where to store data
    fasta_dict = {}
    fasta_dict['ras'] = os.path.join(
        location, '../../data', 'Ras_family_trimmed.fasta'
    )

    return fasta_dict

