from .mutagenesis_visualization import *
import numpy as np
import os
import sys

__author__ = "Frank Hidalgo"
__version__ = "0.0.9"
__title__ = "Mutagenesis Visualization"
__license__ = "GPLv3"
__author_email__ = "fhidalgoruiz@berkeley.edu"


def demo(figure='heatmap'):
    """
    Performs a demonstration of the mutagenesis_visualization software.

    Parameters
    -----------
    figure : str, default 'heatmap'
        There are 5 example plots that can be displayed to test the package is working on your station.
        The 5 options are 'heatmap', 'miniheatmap', 'mean', 'kernel' and 'pca'. Check the documentation for more information.

    Returns
    -------
    None.
    """
    # Use relative file import to access the data folder
    location = os.path.dirname(os.path.realpath(__file__))
    my_file = os.path.join(location, 'data', 'HRas166_RBD.csv')

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')

    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

    # Define secondary structure
    secondary = [['L0'], ['β1']*(9-1), ['L1']*(15-9), ['α1']*(25-15), ['L2']*(36-25), ['β2']*(46-36), ['L3']*(48-46), 
                 ['β3']*(58-48), ['L4'] * (64-58),['α2']*(74-64), ['L5']*(76-74), ['β4']*(83-76), 
                 ['L6']*(86-83), ['α3']*(103-86), ['L7']*(110-103), ['β5']*(116-110), ['L8']*(126-116), ['α4']*(137-126),
                 ['L9']*(140-137), ['β6']*(143-140), ['L10']*(151-143), ['α5']*(172-151), ['L11']*(190-172)]

    # Create object
    hras_RBD = Screen(dataset=hras_enrichment_RBD,
                      sequence=hras_sequence, secondary=secondary)

    if figure == 'heatmap':
        # Create heatmap plot
        hras_RBD.heatmap(title='H-Ras 2-166', show_cartoon=True)
    elif figure == 'miniheatmap':
        # Condensed heatmap
        hras_RBD.miniheatmap(title='Wt residue H-Ras')
    elif figure == 'mean':
        # Mean enrichment by position
        hras_RBD.mean(figsize=[6, 2.5], mode='mean',show_cartoon=True, yscale=[-2, 0.5])
    elif figure == 'kernel':
        # Plot kernel dist using sns.distplot.
        hras_RBD.kernel(histogram=True, title='H-Ras 2-166', xscale=[-2, 1])
    elif figure == 'pca':
        # PCA by amino acid substitution
        hras_RBD.pca(title='', dimensions=[0, 1], figsize=(2, 2), adjustlabels=True)
    return