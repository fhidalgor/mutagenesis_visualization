"""
This module contains the class Screen, which groups the plotting classes.
"""
# Regular libraries
from mutagenesis_visualization.main.heatmaps.heatmap import Heatmap
from typing import List, Optional, Any
import numpy as np
from pandas import DataFrame

from mutagenesis_visualization.main.bar_graphs.mean_barplot import MeanBar
from mutagenesis_visualization.main.bar_graphs.mean_differential import MeanDifferential
from mutagenesis_visualization.main.bar_graphs.mean_position import MeanPosition
from mutagenesis_visualization.main.kernel.kernel import Kernel
from mutagenesis_visualization.main.kernel.histogram import Histogram
from mutagenesis_visualization.main.heatmaps.heatmap import Heatmap
from mutagenesis_visualization.main.heatmaps.heatmap_columns import HeatmapColumns
from mutagenesis_visualization.main.heatmaps.heatmap_rows import HeatmapRows
from mutagenesis_visualization.main.other_stats.box import Box
from mutagenesis_visualization.main.other_stats.cumulative import Cumulative
from mutagenesis_visualization.main.other_stats.miniheatmap import Miniheatmap
from mutagenesis_visualization.main.other_stats.roc_analysis import ROC
from mutagenesis_visualization.main.other_stats.secondary import Secondary

from mutagenesis_visualization.main.pca_analysis.correlation import Correlation
from mutagenesis_visualization.main.pca_analysis.group_correlation import GroupCorrelation
from mutagenesis_visualization.main.pca_analysis.individual_correlation import IndividualCorrelation
from mutagenesis_visualization.main.pca_analysis.pca import PCA
from mutagenesis_visualization.main.plotly.heatmap import HeatmapP
from mutagenesis_visualization.main.plotly.histogram import HistogramP
from mutagenesis_visualization.main.plotly.mean_enrichment import MeanEnrichmentP
from mutagenesis_visualization.main.plotly.rank import RankP
from mutagenesis_visualization.main.plotly.scatter_3d_pdb import Scatter3DPDB
from mutagenesis_visualization.main.plotly.scatter_3d import Scatter3D
from mutagenesis_visualization.main.plotly.scatter import ScatterP
from mutagenesis_visualization.main.scatter.scatter import Scatter
from mutagenesis_visualization.main.utils.snv import select_snv, select_nonsnv
from mutagenesis_visualization.main.utils.pandas_functions import (
    transform_dataset, transform_sequence, transform_secondary
)

try:
    from mutagenesis_visualization.main.pymol.pymol import Pymol
except ModuleNotFoundError:
    pass


class Screen:
    """
    *Screen* represents a mutagenesis experiment. If you are doing deep scan
    mutagenesis, then every amino acid in the protein has been mutated to
    every possible amino acid. For example, if there was a leucine at
    position 2, then this leucine would be mutated to the other 19 naturally
    occurring amino acids. However, you can also use the package if you
    only have a handful of amino acid substitutions.


    Parameters
    -----------
    dataset : array
        2D matrix containing the enrichment scores of the point mutants.
        Columns will contain the amino acid substitutions, rows will
        contain the enrichment for each residue in the protein sequence.

    sequence : str
        Protein sequence in 1 letter code format.

    aminoacids : list, default list('ACDEFGHIKLMNPQRSTVWY*')
        Amino acid substitutions (rows). Submit in the same order that
        is used for the array.

    start_position : int, default 2
        First position in the protein sequence that will be used for the
        first column of the array. If a protein has been mutated only
        from residue 100-150, then if start_position = 100, the algorithm
        will trim the first 99 amino acids in the input sequence. The last
        residue will be calculated based on the length of the input array.
        We have set the default value to 2 because normally the Methionine
        in position 1 is not mutated.

    secondary : list, optional
        This parameter is used to group the data by secondary structure.
        The format is the name of the secondary structure multiplied by
        the residue length of that motif.
        Example : [['β1']*(8),['L1']*(7),['α1']*(9),...,].

    roc_df: Pandas dataframe, optional
        A dataframe that contains a column of variants labeled 'Variant'
        with a column labeled 'Class' containing the true class of that
        mutation. This can be used to compare enrichment scores to some
        label (such as pathogenicity as found in a Cancer database)
        using ROC AUC.

    fillna : float, default 0
        How to replace NaN values.


    Attributes
    ------------
    dataframe : pandas dataframe
        Contains the enrichment scores, position, sequence.

    Other attributes are same as input parameters: dataset, aminoacids,
    start_position, roc_df, secondary

    """
    def __init__(
        self,
        dataset: Any,
        sequence: str,
        aminoacids: List[str] = list('ACDEFGHIKLMNPQRSTVWY*'),
        start_position: Optional[int] = 2,
        fillna: Optional[float] = 0,
        secondary: Optional[list] = None,
        roc_df: Optional[DataFrame] = None
    ):
        # Instances
        self.dataset = np.array(dataset)
        self.aminoacids = aminoacids
        self.start_position = start_position
        self.end_position = len(self.dataset[0]) + start_position
        #self.sequence_raw = ''.join(sequence)  # why I am doing this?
        self.sequence: str = transform_sequence(self.dataset, sequence, self.start_position)
        self.dataframe_stopcodons, self.dataframe = transform_dataset(
            self.dataset, self.sequence, self.aminoacids, self.start_position, fillna
        )
        self.dataframe_snv = select_snv(self.dataframe)
        self.dataframe_nonsnv = select_nonsnv(self.dataframe)

        # Optional parameters
        self.roc_df = roc_df
        self.secondary = secondary
        if self.secondary is not None:
            self.secondary, self.secondary_dup = transform_secondary(
                self.dataset, self.secondary, self.start_position, self.aminoacids
            )

        # Assert messages
        assert len(sequence) >= len(self.dataset[0]), "Input sequence is not long enough."

        # kernel
        self.kernel = Kernel(dataset=self.dataframe['Score_NaN'])
        self.histogram = Histogram(dataframe=self.dataframe, dataframe_snv=self.dataframe_snv, dataframe_nonsnv=self.dataframe_nonsnv)

        # heatmaps
        self.heatmap = Heatmap(
            dataframe = self.dataframe,
            sequence = self.sequence,
            start_position = self.start_position,
            dataframe_stopcodons = self.dataframe_stopcodons,
            secondary = self.secondary
        )

        self.heatmap_rows = HeatmapRows(
            dataframe = self.dataframe,
            sequence = self.sequence,
            start_position= self.start_position,
            dataframe_stopcodons= self.dataframe_stopcodons,
        )
        self.heatmap_columns = HeatmapColumns(
            dataframe= self.dataframe,
            sequence= self.sequence,
            start_position= self.start_position,
            dataframe_stopcodons= self.dataframe_stopcodons,
        )

        # bar
        self.mean = MeanBar(dataframe=self.dataframe)
        #self.differential = MeanDifferential(dataframe=self.dataframe)
        self.position = MeanPosition(dataframe=self.dataframe, start_position = self.start_position)

        # scatter
        #self.scatter = Scatter(self.dataframe)


        # PCA
        self.correlation = Correlation(dataframe_stopcodons = self.dataframe_stopcodons, start_position = self.start_position)
        #self.group_correlation = GroupCorrelation(dataframe_stopcodons = self.dataframe_stopcodons, start_position = self.start_position)
        self.individual_correlation = IndividualCorrelation(dataframe = self.dataframe)
        self.pca = PCA(dataframe = self.dataframe)

        """
        # other stats
        self.box = Box()
        self.cumulative = Cumulative(self.dataframe, self.dataframe_snv, self.dataframe_nonsnv)
        self.miniheatmap = Miniheatmap()
        self.neighboreffect = NeighborEffect()
        self.secondary_mean = Secondary(self.dataframe, self.secondary_dup)
        self.roc = ROC(self.dataframe)
        """
        # plotly
        self.plotly_heatmap = HeatmapP(
            sequence = self.sequence,
            dataframe_stopcodons = self.dataframe_stopcodons)
        self.plotly_histogram = HistogramP(dataframe = self.dataframe)
        self.plotly_mean = MeanEnrichmentP(dataframe = self.dataframe)
        self.plotly_rank = RankP(dataframe = self.dataframe)
        self.plotly_scatter = ScatterP(dataframe = self.dataframe)
        self.plotly_scatter_3D = Scatter3D(dataframe = self.dataframe, start_position = self.start_position,  end_position = self.end_position)
        self.plotly_scatter_3D_pdbprop = Scatter3DPDB(dataframe = self.dataframe, start_position = self.start_position, end_position = self.end_position)

        # pymol
        #try:
            #self.pymol = Pymol(self.dataframe_stopcodons, )
        #except ModuleNotFoundError:
            #pass
