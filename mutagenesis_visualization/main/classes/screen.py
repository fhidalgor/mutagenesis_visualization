"""
This module contains the class Screen, which groups the plotting classes.
"""
# Regular libraries
from typing import List, Optional, Any, Union, Dict
import numpy as np
from pandas import DataFrame
from itertools import combinations
from pathlib import Path

from mutagenesis_visualization.main.bar_graphs.enrichment_bar import EnrichmentBar
from mutagenesis_visualization.main.bar_graphs.differential import Differential
from mutagenesis_visualization.main.bar_graphs.position_bar import PositionBar
from mutagenesis_visualization.main.bar_graphs.secondary import Secondary
from mutagenesis_visualization.main.kernel.kernel import Kernel
from mutagenesis_visualization.main.kernel.histogram import Histogram
from mutagenesis_visualization.main.kernel.multiple_kernels import MultipleKernel
from mutagenesis_visualization.main.heatmaps.heatmap import Heatmap
from mutagenesis_visualization.main.heatmaps.heatmap_columns import HeatmapColumns
from mutagenesis_visualization.main.heatmaps.heatmap_rows import HeatmapRows
from mutagenesis_visualization.main.heatmaps.miniheatmap import Miniheatmap
from mutagenesis_visualization.main.other_stats.rank import Rank
from mutagenesis_visualization.main.other_stats.box import Box
from mutagenesis_visualization.main.other_stats.cumulative import Cumulative
from mutagenesis_visualization.main.other_stats.roc_analysis import ROC
from mutagenesis_visualization.main.pca_analysis.correlation import Correlation
from mutagenesis_visualization.main.pca_analysis.individual_correlation import IndividualCorrelation
from mutagenesis_visualization.main.pca_analysis.pca import PCA
from mutagenesis_visualization.main.plotly.differential import DifferentialP
from mutagenesis_visualization.main.plotly.enrichment_bar import EnrichmentBarP
from mutagenesis_visualization.main.plotly.heatmap import HeatmapP
from mutagenesis_visualization.main.plotly.histogram import HistogramP
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

    fillna : float, default 0
        How to replace NaN values.

    secondary : list, optional
        This parameter is used to group the data by secondary structure.
        The format is the name of the secondary structure multiplied by
        the residue length of that motif.
        Example : [['β1']*(8),['L1']*(7),['α1']*(9),...,].

    replicates: list[np.array, Dataframe], optional
        If you have multiple replicates for that experiment, pass them
        in the same format as dataset.

    Attributes
    ------------
    dataframe : pandas dataframe
        Contains the enrichment scores, position, sequence.

    Other attributes are same as input parameters: dataset, aminoacids,
    start_position, roc_df, secondary

    """
    def __init__(
        self,
        dataset: Union[np.array, DataFrame],
        sequence: str,
        aminoacids: List[str],
        start_position: int = 2,
        fillna: float = 0,
        secondary: Optional[list] = None,
        replicates: Optional[List[Union[np.array, DataFrame]]] = None,
    ):

        if isinstance(dataset, DataFrame):
            self.dataset: np.array = np.array(dataset)
        else:
            self.dataset: np.array = dataset
        if isinstance(aminoacids, list):
            self.aminoacids: List[str] = aminoacids
        elif isinstance(aminoacids, str):
            self.aminoacids: List[str] = list(aminoacids)
        self.start_position: int = start_position
        self.end_position: int = len(self.dataset[0]) + start_position
        self.sequence_raw: str = sequence  # unchanged
        self.sequence: str = transform_sequence(self.dataset, sequence, self.start_position)
        self.dataframe_stopcodons, self.dataframe = transform_dataset(
            self.dataset, self.sequence, self.aminoacids, self.start_position, fillna
        )
        self.dataframe_snv: DataFrame = select_snv(self.dataframe)
        self.dataframe_nonsnv: DataFrame = select_nonsnv(self.dataframe)

        # Optional parameters
        self.secondary: Optional[list] = secondary
        self.secondary_dup: Optional[list] = None
        if self.secondary is not None:
            self.secondary, self.secondary_dup = transform_secondary(
                self.dataset, self.secondary, self.start_position, self.aminoacids
            )

        self.replicates: Optional[List[Union[np.array, DataFrame]]] = replicates
        self.replicates_dataframes: List[DataFrame] = []
        if replicates:
            for replicate in replicates:
                _ , df_replicate = transform_dataset(
                np.array(replicate), self.sequence, self.aminoacids, self.start_position, fillna
            )
            self.replicates_dataframes.append(df_replicate)

        # Assert messages
        assert len(sequence) >= len(self.dataset[0]), "Input sequence is not long enough."

        # kernel
        self.kernel: Kernel = Kernel(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.histogram: Histogram = Histogram(
            dataframe=self.dataframe,
            dataframe_snv=self.dataframe_snv,
            dataframe_nonsnv=self.dataframe_nonsnv, aminoacids=self.aminoacids,
        )
        self.multiple_kernel: MultipleKernel = MultipleKernel(dataframe = self.dataframe, aminoacids=self.aminoacids)

        # heatmaps
        self.heatmap: Heatmap = Heatmap(
            dataframe=self.dataframe,
            sequence=self.sequence,
            start_position=self.start_position,
            dataframe_stopcodons=self.dataframe_stopcodons,
            secondary=self.secondary, aminoacids=self.aminoacids
        )

        self.heatmap_rows: HeatmapRows = HeatmapRows(
            dataframe=self.dataframe,
            sequence=self.sequence,
            start_position=self.start_position,
            dataframe_stopcodons=self.dataframe_stopcodons, aminoacids=self.aminoacids
        )
        self.heatmap_columns: HeatmapColumns = HeatmapColumns(
            dataframe=self.dataframe,
            sequence=self.sequence,
            start_position=self.start_position,
            dataframe_stopcodons=self.dataframe_stopcodons, aminoacids=self.aminoacids
        )
        self.miniheatmap: Miniheatmap = Miniheatmap(
            dataframe=self.dataframe,
            sequence_raw=self.sequence_raw,
            start_position=self.start_position,
            dataframe_stopcodons=self.dataframe_stopcodons,
            dataset=self.dataset,
            aminoacids=self.aminoacids
        )

        # bar
        self.enrichment_bar: EnrichmentBar = EnrichmentBar(
            dataframe=self.dataframe, start_position=self.start_position, aminoacids=self.aminoacids, secondary=self.secondary
        )
        self.differential: Differential = Differential(
            dataframe=self.dataframe, start_position=self.start_position, aminoacids=self.aminoacids, secondary = self.secondary
        )
        self.position_bar: PositionBar = PositionBar(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.secondary_mean: Secondary = Secondary(
            dataframe=self.dataframe, secondary_dup=self.secondary_dup, aminoacids=self.aminoacids
        )

        # scatter
        self.scatter: Scatter = Scatter(dataframe=self.dataframe, aminoacids=self.aminoacids)

        # PCA
        self.correlation: Correlation = Correlation(
            dataframe_stopcodons=self.dataframe_stopcodons, start_position=self.start_position, aminoacids=self.aminoacids
        )
        self.individual_correlation: IndividualCorrelation = IndividualCorrelation(
            dataframe=self.dataframe, aminoacids=self.aminoacids
        )
        self.pca: PCA = PCA(dataframe=self.dataframe, secondary_dup = self.secondary_dup, aminoacids=self.aminoacids)

        # other stats
        self.box: Box = Box(aminoacids=self.aminoacids)
        self.cumulative: Cumulative = Cumulative(
            dataframe=self.dataframe,
            dataframe_snv=self.dataframe_snv,
            dataframe_nonsnv=self.dataframe_nonsnv, aminoacids=self.aminoacids
        )
        self.rank: Rank = Rank(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.roc: ROC = ROC(dataframe=self.dataframe, aminoacids=self.aminoacids)

        # plotly
        self.plotly_differential: DifferentialP = DifferentialP(dataframe=self.dataframe, start_position=self.start_position, aminoacids=self.aminoacids)
        self.plotly_enrichment_bar: EnrichmentBarP = EnrichmentBarP(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.plotly_heatmap: HeatmapP = HeatmapP(
            sequence=self.sequence, dataframe_stopcodons=self.dataframe_stopcodons, aminoacids=self.aminoacids
        )
        self.plotly_histogram: HistogramP = HistogramP(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.plotly_rank: RankP = RankP(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.plotly_scatter: ScatterP = ScatterP(dataframe=self.dataframe, aminoacids=self.aminoacids)
        self.plotly_scatter_3d: Scatter3D = Scatter3D(
            dataframe=self.dataframe,
            start_position=self.start_position,
            end_position=self.end_position, aminoacids=self.aminoacids
        )
        self.plotly_scatter_3d_pdbprop: Scatter3DPDB = Scatter3DPDB(
            dataframe=self.dataframe,
            start_position=self.start_position,
            end_position=self.end_position, aminoacids=self.aminoacids
        )

    def scatter_replicates(self, mode: str = 'pointmutant', output_file: Union[None, str, Path] = None, title: str = None, show: bool=False, close: bool =True) -> None:
        """
        Produce replicates of scatter plots.
        """
        combination_replicates = list(combinations(range(0, len(self.replicates)), 2))
        for comb in combination_replicates:
            x_label: str = "Replicate " + str(comb[0] + 1)
            y_label: str = "Replicate " + str(comb[1] + 1)

            scatter_obj_1: Scatter = Scatter(dataframe=self.replicates_dataframes[comb[0]], aminoacids=self.aminoacids)
            scatter_obj_2: Scatter = Scatter(dataframe=self.replicates_dataframes[comb[1]], aminoacids=self.aminoacids)

            scatter_obj_1(scatter_obj_2,mode = mode, output_file=output_file.with_name(output_file.stem+"_"+str(comb[0])+"_vs_"+str(comb[1])+output_file.suffix), x_label = x_label, y_label = y_label,title = title, show=show, close=close)

        # pymol
        #try:
        #self.pymol: Pymol = Pymol(self.dataframe_stopcodons, )
        #except ModuleNotFoundError:
        #pass
