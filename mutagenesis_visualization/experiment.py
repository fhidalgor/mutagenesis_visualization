from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.plotly.heatmap import HeatmapP
from mutagenesis_visualization.main.kernel.kernel import Kernel
import numpy as np
import pandas as pd
from mutagenesis_visualization.main.heatmaps.heatmap import Heatmap
import numpy as np
from mutagenesis_visualization.main.utils.pandas_functions import (transform_dataset, transform_sequence, transform_secondary)


df = pd.read_csv("mutagenesis_visualization/data/HRas166_RBD.csv")
my_file = "mutagenesis_visualization/data/HRas166_RBD.csv"
# Load enrichment scores
dataset = np.genfromtxt(my_file, delimiter=',')

# Define protein sequence
sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH'

# Define secondary structure
secondary = [['L0'], ['β1'] * (9 - 1), ['L1'] * (15 - 9), ['α1'] * (25 - 15), ['L2'] * (36 - 25),
                ['β2'] * (46 - 36), ['L3'] * (48 - 46), ['β3'] * (58 - 48), ['L4'] * (64 - 58), ['α2'] * (74 - 64),
                ['L5'] * (76 - 74), ['β4'] * (83 - 76), ['L6'] * (86 - 83), ['α3'] * (103 - 86), ['L7'] * (110 - 103),
                ['β5'] * (116 - 110), ['L8'] * (126 - 116), ['α4'] * (137 - 126), ['L9'] * (140 - 137),
                ['β6'] * (143 -  140), ['L10'] * (151 - 143), ['α5'] * (172 - 151), ['L11'] * (190 - 172)]
aminoacids = list('ACDEFGHIKLMNPQRSTVWY*')
start_position = 2
fillna = 0

object = Screen(dataset, sequence, aminoacids, start_position, fillna, secondary)
object.scatter.plot(object)
#object.plotly_scatter_3D_pdbprop.plot("mutagenesis_visualization/data/5p21.pdb")











# Create object
"""sequence = transform_sequence(dataset, sequence, start_position)
dataframe_stopcodons, dataframe = transform_dataset(dataset, sequence, aminoacids, start_position, fillna)
secondary, secondary_dup = transform_secondary(dataset, secondary, start_position, aminoacids)

#heatmap = Heatmap(dataframe, sequence, start_position, dataframe_stopcodons, secondary)
#heatmap.plot(show_cartoon=True)
heatmap_p = HeatmapP(dataframe_stopcodons, sequence)
heatmap_p.plot()"""
